module bimetal

using StrainedBandstructures
using ExtendableFEM
using ExtendableGrids
using SimplexGridFactory
using Triangulate
using GridVisualize
using DrWatson
using DataFrames
using Pardiso
using ExtendableSparse

## start with run_watson() --> result goto data directory
## postprocess data with postprocess(; Plotter = PyPlot) --> images go to plots directory

# configure Watson
@quickactivate "StrainedBandstructures" # <- project name
# set parameters that should be included in filename
watson_accesses = ["scale", "latmis", "femorder", "full_nonlin", "nrefs", "strainm", "mb", "hz", "grid_type", "bc", "mstruct","stressor_x"]
watson_allowedtypes = (Real, String, Symbol, Array, DataType)
watson_datasubdir = "bimetal_watson"


function run_watson(; dim = 2, force::Bool = false, generate_vtk = true)

    allparams = Dict(
        "latmis" => Array{Float64,1}(0:0.02:0.2),                                           # lattice mismatch between lattice constants of material A and B (lc = [5,5*(1+latmis)])
        "E" => [[1e-6, 1e-6]],                                                              # elastic moduli of material A and B
        "ν" => [[0.15, 0.15]],                                                              # Poisson numbers of material A and B
        "strainm" => [dim == 2 ? NonlinearStrain2D : NonlinearStrain3D],                    # strain model
        "full_nonlin" => [true],                                                            # use complicated model (ignored if linear strain is used)
        "use_emb" => [true],                                                                # use embedding (true) or damping (false) solver ?
        "nsteps" => [4],                                                                    # number of embedding steps in embedding solver
        "maxits" => [100],                                                                  # max number of iteration in each embedding step
        "tres" => [1e-12],                                                                  # target residual in each embedding step
        "scale" => dim == 2 ? [[50,2000], [100, 2000]] : [[50,50,2000], [100,100,2000]],    # dimensions of bimetal
        "mb" => [0.5],                                                                      # share of material A vs. material B
        "femorder" => [2],                                                                  # order of the finite element discretisation
        "upscaling" => [0],                                                                 # upscaling of results (does so many extra nrefs for plots)
        "nrefs" => [1],                                                                     # number of uniform refinements before solve
        "avgc" => [2],                                                                      # lattice number calculation method (average case)
        "linsolver" => ExtendableSparse.MKLPardisoLU,                                       # linear solver (try ExtendableSparse.MKLPardisoLU or ExtendableSparse.LUFactorization)
        "grid_type" => "default",                                                           # grid options: default, condensator, condensator_tensorgrid
        "bc" => "robin",                                                                    # boundary conditions: robin (Dirichlet and/or Newmann), periodic
    )

    dicts = dict_list(allparams)
    @info "Starting bimetal simulations..." dicts

    # run the simulations and save data
    for (i, d) in enumerate(dicts)
        run_single(d; force = force, generate_vtk = generate_vtk)
    end
end


# defaults for use with run_single
function get_defaults()
    params = Dict(
        "latmis" => 0.05,                               # lattice mismatch between lattice constants of material A and B (lc = [5,5*(1+latmis)])
        "E" => [1e-6, 1e-6],                            # elastic moduli of material A and B
        "ν" => [0.15, 0.15],                            # Poisson numbers of material A and B
        "scale" => [50,50,2000],                        # dimensions of bimetal
        "full_nonlin" => true,                          # use complicated model (ignored if linear strain is used)
        "use_emb" => true,                              # use embedding (true) or damping (false) solver ?
        "nsteps" => 4,                                  # number of embedding steps in embedding solver
        "maxits" => 20,                                 # max number of iteration in each embedding step
        "tres" => 1e-12,                                # target residual in each embedding step
        "mb" => 0.5,                                    # share of material A vs. material B
        "hz" => 50,                                     # lenght in nn between z cuts
        "femorder" => 2,                                # order of the finite element discretisation
        "upscaling" => 0,                               # upscaling of results (does so many extra nrefs for plots)
        "nrefs" => 1,                                   # number of uniform refinements before solve
        "avgc" => 2,                                    # lattice number calculation method (average case)
        "linsolver" => ExtendableSparse.MKLPardisoLU,   # linear solver (try ExtendableSparse.MKLPardisoLU or ExtendableSparse.LUFactorization)
        "grid_type" => "default",                       # grid options: default, condensator, condensator_tensorgrid
        "bc" => "robin",                                # boundary conditions: robin (Dirichlet and/or Newmann), periodic
        "scenario" => 1,                                # scenario number that fixes materials for core/stressor (1: default, no specific material, 2: semiconductor materials)
        "stressor_x" => 0.5,                            # x value for x-dependent stressor material
        "mstruct" => ZincBlende001,                     # material structure type
        "strainm" => NonlinearStrain3D,                 # strain model
        "estrainm" => IsotropicPrestrain,               # elastic strain model
    )
    return params
end
function set_params!(d; kwargs)
    for (k,v) in kwargs
        d[String(k)]=v
    end
    d["strainm"] = length(d["scale"]) == 2 ? NonlinearStrain2D : NonlinearStrain3D
    return nothing
end

function get_data(d = nothing; kwargs)
    ## load parameter set
    if d === nothing
        d = get_defaults()
    end
    set_params!(d; kwargs)
    return d
end


function run_single(d = nothing; force::Bool = false, generate_vtk = true, Plotter = nothing, kwargs...)

    d = get_data(d; kwargs)
    @show d

    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    if isfile(datadir(watson_datasubdir, filename)) && !force
        @info "Skipping dataset $filename... (run with force = true to enforce recomputation)"
        return nothing
    else
        @info "Running dataset $filename..."
    end

    fulld = copy(d)

    if d["scenario"] == 1
        ## compute lattice_mismatch
        lc = [5, 5 * ( 1 + d["latmis"] )]
        @info "lattice mismatch: $(round(100*(lc[2]/lc[1]-1),digits=3)) %"
    elseif d["scenario"] == 2
        materials = [GaAs, AlInAs{d["stressor_x"]}]
        #materials = [GaN, InGaN{d["stressor_x"]}]
        materialstructuretype = d["mstruct"]
        MD = HeteroStructureData(materials, materialstructuretype)

        ## compute lattice_mismatch
        if d["latmis"] == nothing
            nregions = length(MD.data)
            lc::Array{Float64,1} = zeros(Float64,nregions)
            for j = 1 : nregions
                lc[j] = MD.data[j].LatticeConstants[1]
            end
        else
            lc = [5, 5 * ( 1 + d["latmis"] )]
        end
        @info "lattice mismatch: $(round(100*(lc[2]/lc[1]-1),digits=3)) %"
        @info lc

        ## save data
        fulld["MD"] = MD
    else
        @error "scenario not defined"
    end

    ## compute misfit strain
    mb = d["mb"]
    avgc = d["avgc"]
    scale = d["scale"]
    misfit_strain, fulld["α"] = get_lattice_mismatch_bimetal(avgc, [scale[1] * mb, scale[1] * (1 - mb)], lc)
    # fulld["misfit_strain"] = [[misfit_strain[1],misfit_strain[1],-2*misfit_strain[1]*MD.data[1].ElasticConstants["C12"]/MD.data[1].ElasticConstants["C11"]],
    #                            [misfit_strain[2],misfit_strain[2],-2*misfit_strain[2]*MD.data[2].ElasticConstants["C12"]/MD.data[2].ElasticConstants["C11"]]]
    fulld["misfit_strain"] = [misfit_strain[1]*ones(Float64,3), misfit_strain[2]*ones(Float64,3)]
    #fulld["misfit_strain"] = [[0.,0.,0],[0.021986187362812704,0.021986187362812704,0.021986187362812704]]
    @info fulld["misfit_strain"]
    @info fulld["α"]

    ## run simulation
    solution, residual = main(fulld; Plotter = Plotter)

    ## save additional data
    fulld["solution"] = solution
    fulld["residual"] = residual

    ## save data to file
    wsave(datadir(watson_datasubdir, filename), fulld)

    if generate_vtk
        export_vtk(fulld; upscaling = fulld["upscaling"], kwargs...)
    end

    return fulld
end

## this calculates with user-given misfit strain
function main(d::Dict; Plotter = Plotter, verbosity = 0)

    ## unpack paramers
    @unpack scenario, linsolver, latmis, misfit_strain, α, full_nonlin, use_emb, nsteps, maxits, tres, scale, mb, hz, femorder, nrefs, strainm, estrainm, avgc, grid_type, bc = d
    
    if scenario == 1
        ## compute Lame' coefficients μ and λ from ν and E
        @unpack E, ν = d
        μ = E ./ (2  .* (1 .+ ν))
        λ = E .* ν ./ ( (1 .- 2*ν) .* (1 .+ ν))
        nregions = 2
    elseif scenario == 2
        @unpack MD = d
        nregions = length(MD.data)
        C11::Array{Float64,1} = zeros(Float64,nregions)
        C12::Array{Float64,1} = zeros(Float64,nregions)
        for j = 1 : nregions
            C11[j] = MD.data[j].ElasticConstants["C11"]
            C12[j] = MD.data[j].ElasticConstants["C12"]
        end

        ## compute Lame' coefficients μ and λ from elasticity matrix C and update E and ν
        λ = 1e-9 .* C12
        μ = 1e-9 .* (C11 .- C12) ./ 2

        @info λ
        @info μ

        d["E"] = μ .* (3 .* λ .+ 2 .* μ) ./ (μ .+ λ)
        d["ν"] = λ ./ (2 .* (λ .+ μ))
    else
        @error "scenario not defined"
    end
    parameters::Array{Float64,1} = [full_nonlin ? 1 : 0]

    ## generate bimetal mesh
    dim::Int = length(scale)
    periodic_regions = []
    @assert (femorder in 1:3)
     ### DOUBLE CHECK dirichlet conditions ###
    if dim == 3
        if grid_type == "default"
            #xgrid = bimetal_strip3D(; material_border = mb, scale = scale)
            #xgrid = uniform_refine(xgrid,nrefs)
            #xgrid = bimetal_tensorgrid(; scale = scale, nrefs = nrefs, material_border = mb); dirichlet_regions = [3,4]

            xgrid, xgrid_cross_section = bimetal_tensorgrid_uniform(; scale = scale, nrefs = nrefs, material_border = mb, hz = hz); dirichlet_regions = [7,8]

        elseif grid_type == "condensator"
            xgrid = condensator3D(; scale = scale, d = 10, nrefs = nrefs); dirichlet_regions = [1,2,5,6] # core sides and bottoms
        elseif grid_type == "condensator_tensorgrid"
            xgrid, xgrid_cross_section = condensator3D_tensorgrid!(; scale=scale, d=3, dx=1, nrefs=nrefs, stressor_cell_per=10)
            # core sides = [1,2,3,4], bottom and upper sides = [5,6], stressor sides = [7,8,9,10]

            #xgrid, xgrid_cross_section = condensator3D_tensorgrid(; scale = scale, d = 10, nrefs = nrefs)
            dirichlet_regions = [] # free boundary (up to normal direction on bottom + 1 fixed point)
            #dirichlet_regions = [5,6] # bottom and top
            #dirichlet_regions = [7] # stressor side
            #dirichlet_regions = [1,2,3,4,5,6] # core sides and bottoms
            if bc =="periodic"
                periodic_regions = [[1,3,(f1,f2) -> abs(f1[1] - f2[1]) + abs(f1[3] - f2[3]) < 1e-12,  [-1,-1,-1]],
                                [2,4,(f1,f2) -> abs(f1[2] - f2[2]) + abs(f1[3] - f2[3]) < 1e-12,  [-1,-1,-1]],
                                [7,9,(f1,f2) -> abs(f1[1] - f2[1]) + abs(f1[3] - f2[3]) < 1e-12,  [-1,-1,-1]],
                                [8,10,(f1,f2) -> abs(f1[2] - f2[2]) + abs(f1[3] - f2[3]) < 1e-12,  [-1,-1,-1]]]
            end
         else
            xgrid = bimetal_strip3D_middle_layer(; scale = scale, reflevel = nrefs); dirichlet_regions = [3,4]
        end
    else
        if grid_type == "default"
            xgrid = bimetal_strip2D(; material_border = mb, scale = scale)
            dirichlet_regions = [1]
        elseif grid_type == "condensator"
            xgrid = condensator2D_periodic(; A = scale[1], B = scale[2], d = 5, reflevel = nrefs)
            dirichlet_regions = [1,3]
        end
        if bc =="periodic"
            periodic_regions = [[5,6,(f1,f2) -> abs(f1[1] - f2[1]) + abs(f1[2] - f2[2]) < 1e-12,  [-1,-1]],
                            [1,3,(f1,f2) -> abs(f1[1] - f2[1]) + abs(f1[3] - f2[3]) < 1e-12,  [-1,-1]],
                            [2,4,(f1,f2) -> abs(f1[3] - f2[3]) + abs(f1[2] - f2[2]) < 1e-12,  [-1,-1]]]
        end
        xgrid = uniform_refine(xgrid,nrefs)
    end
    if Plotter !== nothing
        gridplot(xgrid_cross_section; Plotter=Plotter)
        @show xgrid
    end

    ## setup model
    if strainm <: NonlinearStrain3D && dim == 2
        @warn "changing strainm to $strainm due to non-matching dimension"
        strainm = NonlinearStrain2D
        d["strainm"] = strainm
    elseif strainm <: LinearStrain3D && dim == 2
        @warn "changing strainm to $strainm due to non-matching dimension"
        strainm = LinearStrain2D
        d["strainm"] = strainm
    end
    full_nonlin *= strainm <: NonlinearStrain

    ###########################
    ### PROBLEM DESCRIPTION ###
    ###########################

    ## create PDEDescription and add displacement unknown
    Problem = ProblemDescription("nanowire bending")
    u = Unknown("u"; name = "displacement")
    assign_unknown!(Problem, u)

    ## add (nonlinear) operators for displacement equation
    if scenario == 1
        DisplacementOperator = get_displacement_operator([IsotropicElasticityTensor(λ[r], μ[r], dim) for r = 1 : nregions], strainm, IsotropicPrestrain, misfit_strain, α; dim = dim, displacement = u, emb = parameters, regions = 1:nregions, bonus_quadorder = 2*(femorder-1))
    elseif scenario == 2
        DisplacementOperator = get_displacement_operator(MD.TensorC, strainm, IsotropicPrestrain, misfit_strain, α; dim = dim, displacement = u, emb = parameters, regions = 1:nregions, bonus_quadorder = bonus_quadorder)
    end
    assign_operator!(Problem, DisplacementOperator)

    ## choose finite element type for displacement (vector-valued) and polarisation (scalar-valued)
    if femorder == 1
        FEType_D = H1P1{dim}
        FEType_P = H1P1{1}
    elseif femorder == 2
        FEType_D = H1P2{dim,dim}
        FEType_P = H1P2{1,dim}
    elseif femorder == 3
        FEType_D = H1P3{dim,dim}
        FEType_P = H1P3{1,dim}
    end

    ## add Dirichlet boundary data on front
    BoundaryOperator = nothing
    PeriodicBoundaryOperator = nothing
    if bc != "periodic"
        if length(dirichlet_regions) > 0
            BoundaryOperator = HomogeneousBoundaryData(u; regions = dirichlet_regions)
            assign_operator!(Problem, BoundaryOperator)
        end
    elseif bc == "periodic"
        ## add coupling information for periodic regions
        FES = FESpace{FEType_D}(xgrid)
        dofsX, dofsY, factors = Int[], Int[], Int[]
        for j = 1 : length(periodic_regions)
            dofsX_j, dofsY_j, factors_j = get_periodic_coupling_info(FES, xgrid, periodic_regions[j][1], periodic_regions[j][2], periodic_regions[j][3]; factor_components = periodic_regions[j][4])
            append!(dofsX, dofsX_j)
            append!(dofsY, dofsY_j)
            append!(factors, factors_j)
        end
        PeriodicBoundaryOperator = CombineDofs(1, 1, dofsX, dofsY, factors)
        assign_operator!(Problem, PeriodicBoundaryOperator)
   end

    ##############
    ### SOLVER ###
    ##############

    polarisation = false
    PolarisationOperator = nothing

    Solution, residual = solve_lowlevel(xgrid,
                                BoundaryOperator,
                                DisplacementOperator.kernel,
                                PolarisationOperator,
                                parameters;
                                periodic_boundary_operator = PeriodicBoundaryOperator,
                                linsolver = linsolver,
                                nsteps = [nsteps, 1],
                                FETypes = [FEType_D, FEType_P],
                                target_residual = [tres, tres],
                                solve_polarisation = polarisation,
                                coupled = false,
                                maxiterations = [maxits, 1])

    return Solution, residual
end

function get_lattice_mismatch_bimetal(avgc, geometry, lc)
    r::Array{Float64,1} = geometry[1:2]
    a::Array{Float64,1} = zeros(Float64,2)

    A_core = r[1]
    A_stressor = r[2]
    lc_avg = (lc[1]*A_core + lc[2]*A_stressor)/(A_core + A_stressor)

    for j = 1 : 2
        if avgc == 1
            a[j] = (lc[j] - lc[1])/lc[j]
        elseif avgc == 2
            a[j] = (lc[j] - lc_avg)/lc[j]
        elseif avgc == 3
            a[j] = (lc[j] - lc_avg)/lc_avg
        end
    end

    #return a .* (1 .+ a./2), a
    return a, a
end


## load a result dictionary (to browse results in julia)
function load_data(d = nothing; kwargs)
    ## complete parameter set
    d = get_data(d; kwargs)

    # load dict from file
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    @info filename
    d = wload(datadir(watson_datasubdir, filename))
    return d
end

function export_vtk(d = nothing; upscaling = 0, kwargs...)
    d = load_data(d; kwargs)
    filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    @unpack solution, strainm, estrainm, misfit_strain = d
    repair_grid!(solution[1].FES.xgrid)
    StrainedBandstructures.exportVTK(datadir(watson_datasubdir, filename_vtk), misfit_strain, solution[1], nothing; EST = estrainm, strain_model = strainm, P0strain = true, upscaling = upscaling)
end

function export_cuts(; 
    scale = [50, 50, 2000],
    nrefs = 1,
    upscaling = 0,
    eps_gfind = 1e-11,
    mb = 0.5,
    cut_direction = 1, # 1 = x-direction, 2 = y, 3 = z
    femorder = "max",
    strainm = length(scale) == 3 ? NonlinearStrain3D : NonlinearStrain2D,
    cut_levels = [scale[3] * 0.5],
    Plotter = nothing)

    @info "Exporting cuts..."

    # load all data
    alldata = collect_results(datadir(watson_datasubdir);subfolders=true)

    # filter and sort
    df = filter(:scale => ==(scale), alldata)
    df = filter(:nrefs => ==(nrefs), df)
    df = filter(:mb => ==(mb), df)
    df = filter(:strainm => ==(strainm), df)

    if femorder == "max"
        femorder = maximum(df[!,:femorder])
    end
    if upscaling == "auto"
        upscaling::Int = femorder > 1
    end

    df = filter(:femorder => ==(femorder), df)
    @show df

    # compute cuts
    for data in eachrow(df)
        solution = data[:solution]
        strainm = data[:strainm]
        eps0 = data[:misfit_strain]

        repair_grid!(solution[1].FES.xgrid)
        xgrid_plot = solution[1].FES.xgrid
        
        Solution_plot = nothing
        if upscaling > 0
            xgrid_plot = uniform_refine(xgrid_plot, upscaling; store_parents = true)
            if length(solution) == 1
                FES = FESpace{eltype(solution[1].FES)}(xgrid_plot)
                Solution_plot = FEVector{Float64}("u_h (upscale)", FES)
            else
                FES = [FESpace{eltype(solution[1].FES)}(xgrid_plot), FESpace{eltype(solution[2].FES)}(xgrid_plot)]
                Solution_plot = FEVector{Float64}(["u_h (upscale)", "V_P (upscale)"],[FES[1], FES[2]])
            end
            for j = 1 : length(solution)
                interpolate!(Solution_plot[j], solution[j]; use_cellparents = true, eps = eps_gfind)
            end
        else
            Solution_plot = solution
        end

        filename_cuts = savename(data, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses) * "_CUTS/"
        mkpath(datadir(watson_datasubdir, filename_cuts))
        plane_points = [[0.25*scale[1],0.25*scale[2]],[0.75*scale[1],0.25*scale[2]],[0.25*scale[1],0.75*scale[2]]] # = [X,Y] coordinates of the three points that define the cut plane
        perform_simple_plane_cuts(datadir(watson_datasubdir, filename_cuts), Solution_plot, plane_points, cut_levels; eps0 = eps0, cut_direction = cut_direction, eps_gfind = 1e-10, only_localsearch = true, strain_model = strainm, Plotter = Plotter, export_uniform_data = true)
    end
end

function postprocess(; 
    scales = [[50, 2000], [100, 2000]],
    maxlc = [0.1, 0.2],
    nrefs = 1,
    hz = 50,
    mb = 0.5,
    femorder = "max",
    strainm = length(scales[1]) == 3 ? NonlinearStrain3D : NonlinearStrain2D,
    Plotter = nothing)

    @assert Plotter !== nothing "need a Plotter (e.g. PyPlot)"
    Plotter.close("all")
    @info "Starting postprocessing..."

    # load all data
    alldata = collect_results(datadir(watson_datasubdir);subfolders=true)

    # init plot
    fig, (ax1, ax2) = Plotter.subplots(2, 1);
    marker = ["o","s"]
    color = ["red", "blue"]
    legend = []

    for j = 1 : length(scales)

        # filter and sort
        scale = scales[j]
        df = filter(:scale => ==(scale), alldata)
        df = filter(:nrefs => ==(nrefs), df)
        df = filter(:hz => ==(hz), df)
        df = filter(:mb => ==(mb), df)
        df = filter(:strainm => ==(strainm), df)
        df = filter(:latmis => <=(maxlc[j]), df)

        if femorder == "max"
            femorder = maximum(df[!,:femorder])
        end
        df = filter(:femorder => ==(femorder), df)
        sort!(df, [:latmis])
        @show df

        # compute curvature etc
        lattice_mismatch = []
        sim_angle = []
        ana_angle = []
        sim_curvature = []
        ana_curvature = []
        for data in eachrow(df)
            solution = data[:solution]
            scale = data[:scale]
            # misfit_strain = data[:misfit_strain]
            ## compute bending statistics
            scaling = 1 - abs(2*mb-1)
            bending_axis_end_points = [[scaling*scale[1],scale[2]/2,0],[scaling*scale[1],scale[2]/2,scale[3]]]
            @info bending_axis_end_points
            angle, curvature, dist_bend, farthest_point = compute_statistics(solution[1].FES.xgrid, solution[1], bending_axis_end_points, eltype(solution[1].FES))

            # calculate analytic curvature (is it correct for 3D?)
            E = data[:E]
            α = data[:α]
            mb = data[:mb]
            #factor = 1/2*(α[2] - α[1])*(2 + α[1] + α[2])
            factor = α[2] - α[1]
            n = E[1]/E[2]
            h1 = data[:scale][1] * (1-mb)
            h2 = data[:scale][1] * (mb)
            m = h1/h2
            analytic_curvature = abs( 6.0 * factor * (1+m)^2 / ((h1+h2) * ( 3*(1+m)^2 + (1+m*n)*(m^2+1/(m*n)))) )
            analytic_arc = scale[3]*(m^4 + 4*m*n + 6*m^2*n + 4*m^3*n + n^2 + m*(m^3 + 4*n + 3*m*n)*α[1] + n*(3*m^2 + 4*m^3 + n)*α[2])/(m^4 + 4*m*n + 6*m^2*n + 4*m^3*n + n^2)
            analytic_arc = mod(analytic_arc,2*π/analytic_curvature)
            analytic_angle = mod(analytic_arc*analytic_curvature/2, 2*π) * 180/π
            # if farthest_point[3] > 0
            #     analytic_angle = asin(dist_bend/2*analytic_curvature) * 180/π
            # else
            #     analytic_angle = 180 - asin(dist_bend/2*curvature) * 180/π
            # end
            #@info "dist_bend = $(dist_bend)"
            @info "simulation ===> R = $(1/curvature) | curvature = $curvature | bending angle = $(angle)°"
            @info "analytic   ===> R = $(1/analytic_curvature) | curvature = $analytic_curvature | bending angle = $(analytic_angle)°"
            println()
            push!(lattice_mismatch, data[:latmis])
            push!(sim_angle, angle)
            push!(ana_angle, analytic_angle)
            push!(sim_curvature, curvature)
            push!(ana_curvature, analytic_curvature)
        end
        for i in eachindex(sim_curvature)
            println(sim_curvature[i])
        end

        @info sim_curvature

        @info "Plotting ..."
        lattice_mismatch *= 100
        ax1.plot(lattice_mismatch, 1e3 .* sim_curvature, color = color[j], marker = marker[j])
        ax1.plot(lattice_mismatch, 1e3 .* ana_curvature, color = color[j], linestyle = "--")
        ax2.plot(lattice_mismatch, sim_angle, color = color[j], marker = marker[j])
        ax2.plot(lattice_mismatch, ana_angle, color = color[j], linestyle = "--")

        ax1.set_title("Stressor: $(Int(mb*100))%, Core: $(Int((1-mb)*100))%")
        append!(legend, ["simulation, d = $(scale[1])","analytic, d = $(scale[1])"])
        ax1.set_ylabel("curvature (μm^-1)")
        ax2.set_ylabel("angle (degrees)")
        ax2.set_xlabel("lattice mismatch (%)")
    end

    ax1.legend(legend)
    ax2.legend(legend)
    Plotter.savefig("plots/curvature_angle_material_border=$(mb).png")
end

end