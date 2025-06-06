module nanowire

using StrainedBandstructures
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using Printf
using SimplexGridFactory
using Triangulate
using LinearAlgebra
using TetGen
using DrWatson
using DataFrames
using LinearSolve
using Metis
using Pardiso

## start with run_watson() --> result goto data directory
## postprocess data with postprocess(; Plotter = PyPlot) --> images go to plots directory

# configure Watson
@quickactivate "StrainedBandstructures" # <- project name
# set parameters that should be included in filename
watson_accesses = ["geometry", "nrefs", "femorder", "mstruct", "stressor_x", "full_nonlin", "interface_refinement", "cross_section", "rotate"]
watson_allowedtypes = (Real, String, Symbol, Array, DataType)
watson_datasubdir = "nanowire"


## get default parameters
function get_defaults()
    params = Dict(
        "shell_x" => 0.3,                               # x value for x-dependent shell material
        "stressor_x" => 0.5,                            # x value for x-dependent stressor material
        "strainm" => NonlinearStrain3D,                 # strain model
        "estrainm" => IsotropicPrestrain,               # elastic strain model
        "full_nonlin" => true,                          # use complicated model (ignored if linear strain is used)
        "use_emb" => true,                              # use embedding (true) or damping (false) solver ?
        "nsteps" => 5,                                  # number of embedding steps in embedding solver
        "maxits" => 20,                                 # max number of iteration in each embedding step
        "tres" => 1e-12,                                # target residual in each embedding step
        "geometry" => [30, 20, 15, 2000],               # dimensions of nanowire
        "scenario" => 1,                                # scenario number that fixes materials for core/shell/stressor
        "mb" => 0.5,                                    # share of material A vs. material B
        "mstruct" => ZincBlende001,                     # material structure type
        "femorder" => 1,                                # order of the finite element discretisation (displacement)
        "use_lowlevel_solver" => true,                  # use new implementation based on low level structures (should be faster)
        "bonus_quadorder" => 2,                         # bonus quadrature order for nonlinear operator assembly
        "femorder_P" => 1,                              # order of the finite element discretisation (polarisation)
        "nrefs" => 0,                                   # number of uniform refinements before solve
        "avgc" => 1,                                    # lattice number calculation method (average case)
        "polarisation" => false,                         # also solve for polarisation
        "fully_coupled" => false,                       # parameter for later when we have the full model
        "postprocess" => false,                         # angle calculation, vtk files, cuts
        "linsolver" => LinearSolve.MKLPardisoFactorize(),         # linear solver (everything supported by LinearSolve)
        "npartitions" => 1,                             # number of partitions for parallel solve and assembly
        "parallel" => false,                             # number of partitions for parallel solve and assembly
        "damping" => 0,                                 # damping in Newton iteration
        "uniform_grid" => true,                         # uniform grid in z direction or nonuniform grid with local refinement at cut levels
        "interface_refinement" => false,                # enables a finer grid at material interface for each cross-section
        "z_levels_dist" => 100,                         # distance between z-levels is z_levels_dist * 2^(-nrefs)
        "cut_levels_pos" => [0.5],                      # position of cut-levels w.r.t. nanowire length. i.e., cut_levels = cut_levels_pos * geometry[4]
        "cross_section" => 1,                           # geometry of stressor; either version 1,2 or 3
        "corner_refinement" => false,                   # assigning more nodes at interface corners
        "rotate" =>  0                                  # rotate cross section clockwise
    )
    return params
end

function set_params!(d; kwargs...)
    for (k,v) in kwargs
        d[String(k)]=v
    end
    return nothing
end


function get_scenario(scenario, shell_x, stressor_x, materialstructuretype)
    if scenario == 1
        materials = [GaAs,GaAs,AlInAs{stressor_x}]                  # scenario 1
    elseif scenario == 2
        materials = [GaAs,AlInAs{shell_x},AlInAs{stressor_x}]       # scenario 2
    elseif scenario == 3
        materials = [GaAs,GaAs,AlGaAs{stressor_x}]                  # scenario 3
    elseif scenario == 4
        materials = [GaAs,GaAs,InGaAs{stressor_x}]                  # scenario 4
    elseif scenario == 5
        materials = [InGaAs{0.53},InGaAs{0.53},AlInAs{stressor_x}]  # scenario 5
    else
        @error "scenario not defined"
    end
    MD = HeteroStructureData(materials, materialstructuretype)
    return MD
end

# call this function for the full nanowire scenario
# default parameters (see above) can be modified by kwargs
# Plotter is used for postprocessing
function main(d = nothing; verbosity = 0, Plotter = nothing, force::Bool = false, old_grid = false, kwargs...)

    ## load parameter set
    if d === nothing
        d = get_defaults()
        set_params!(d; kwargs...)
    end

    ## raise a warning if anisotropic material is used with isotropic elastic strain_model
    if d["mstruct"] != ZincBlende001 && d["estrainm"] == IsotropicPrestrain
        @error "ArgumentError: Use IsotropicPrestrain only with ZincBlende001 materials"
        return
    end

    println("***Solving nanowire problem***")

    ## set log level
    #set_verbosity(verbosity)

    ################
    ### SCENARIO ###
    ################
    ## choose material data DofMapTypes
    @unpack shell_x, stressor_x, scenario, mstruct = d
    MD = get_scenario(scenario, shell_x, stressor_x, mstruct)
    nregions = length(MD.data)
    println("SCENARIO DATA")
    println("=============")
    for j = 1 : length(MD.data)
        println("        material[$j] = $(String(get_materialtype(MD.data[j])))")
    end
    println("          structure = $mstruct")

    #######################
    ### MODEL PARAMETER ###
    #######################
    @unpack geometry, fully_coupled, full_nonlin, strainm, estrainm, avgc, femorder, femorder_P, polarisation = d
    @assert fully_coupled == false "fully coupled model not yet implemented"
    ## setup parameter Array (which is also returned and can be used later to embedd)
    parameters::Array{Float64,1} = [full_nonlin ? 1 : 0]
    ## compute lattice misfit strain
    eps0, a = get_lattice_misfit_nanowire(avgc, MD, geometry, full_nonlin)
    @show eps0
    @show a

    println("        strain type = $(strainm)")
    println("   elasticity model = $(full_nonlin ? "full" : "simplified")")
    println("      misfit strain = $(eps0)")
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    println("           filename = $filename\n")

    # check if result files already exist
    if isfile(datadir(watson_datasubdir, filename)) && !force
        @warn "Result files already exist, rerun with force = true if needed"
        return nothing
    end

    ############
    ### GRID ###
    ############
    @unpack nrefs = d
    ## load mesh, expecting the following regions
    ##              regions [1,2,3] = [core, shell, stressor]
    ##      boundaryregions [1,2,3] = [core front, shell boundary, stressor boundary]
    if old_grid
        gridfile = "grids/nanowire-grid($(geometry[1]),$(geometry[2]),$(geometry[3]),$(geometry[4])).sg" # default: gridfile = "nanowire-grid(30,20,15,2000).sg"
        xgrid = simplexgrid(gridfile)
        xgrid = uniform_refine(xgrid,nrefs)
    else
        geometry = Array{Float64,1}(geometry)
        if d["cross_section"] == 1 || d["cross_section"] == 2
            geometry[1] /= sqrt(3)  # core hexagon side
            geometry[2] /= sqrt(3)  # shell hexagon side
                                    # geometry[3] is the stressor width
        elseif d["cross_section"] == 3
            # geometry[1] and geometry[2] are the core/shell hexagon sides, i.e., geometry[1]+geometry[2] is the diameter of the core+shell hexagon
            geometry[3] *= 2/sqrt(3) # geometry[3] is the stressor hexagon side
        end
        if d["uniform_grid"] == true
            d["cut_levels"] = nothing
        else
            d["cut_levels"] = d["cut_levels_pos"]*d["geometry"][4]   # z-level locations of cuts

        end
        if d["interface_refinement"] == true
            if d["cross_section"] == 1
                refinement_width = 1/4*geometry[3]
            elseif d["cross_section"] == 2
               refinement_width = (geometry[3] >= 10 ? 4 : 1 + (geometry[1]+geometry[2])/25)/sqrt(3)
            elseif d["cross_section"] == 3
                refinement_width = 1/4*geometry[2]
            end
           manual_refinement = true
        else
            refinement_width = nothing
            manual_refinement = false
        end
        xgrid, xgrid_cross_section = nanowire_tensorgrid_mirror(; scale = geometry, nrefs = nrefs, z_nrefs = 2, shape = d["cross_section"],
                                            cut_levels = d["cut_levels"], z_levels_dist = d["z_levels_dist"],
                                            refinement_width = refinement_width, corner_refinement = d["corner_refinement"],
                                            manual_refinement = manual_refinement, rotate = d["rotate"])
    end

    npartitions = d["npartitions"]
    if npartitions > 1
        xgrid = partition(xgrid, PlainMetisPartitioning(npart = npartitions))
    end

    #xgrid = nanowire_grid(; scale = geometry)
    gridplot(xgrid_cross_section; Plotter=Plotter)
    @show xgrid


    ##############
    ### SOLVER ###
    ##############
    regions_bc = [4,5,6]  # 4 = core bottom, 5 = shell bottom, 6 = stressor bottom

    ## choose finite element type for displacement (vector-valued) and polarisation (scalar-valued)
    if femorder == 1
        FEType_D = H1P1{3}
        FEType_P = H1P1{1}
    elseif femorder == 2
        FEType_D = H1P2{3,3}
        FEType_P = H1P2{1,3}
    elseif femorder == 3
        FEType_D = H1P3{3,3}
        FEType_P = H1P3{1,3}
    end

    ## call solver
    @unpack nsteps, tres, maxits, linsolver, bonus_quadorder, damping, use_lowlevel_solver = d

    ## define unknnowns
    u = Unknown("u"; name = "discplacement")
    V = Unknown("V"; name = "Polarisation potential")

    if (use_lowlevel_solver)
        @assert polarisation == false
        BoundaryOperator = HomogeneousBoundaryData(1; regions = regions_bc)
        DisplacementOperator = get_displacement_operator(MD.TensorC, strainm, estrainm, eps0, a; dim = 3, displacement = u, emb = parameters, regions = 1:nregions, bonus_quadorder = bonus_quadorder)
        PolarisationOperator = nothing
        if d["parallel"]
            Solution, last_residual = solve_lowlevel_parallel(xgrid,
                                    BoundaryOperator,
                                    DisplacementOperator.kernel,
                                    PolarisationOperator,
                                    parameters;
                                    linsolver = linsolver,
                                    nsteps = [nsteps, 1],
                                    FETypes = [FEType_D, FEType_P],
                                    target_residual = [tres, tres],
                                    solve_polarisation = polarisation,
                                    coupled = false,
                                    maxiterations = [maxits, 1])
        else
            Solution, last_residual = solve_lowlevel(xgrid,
                                    BoundaryOperator,
                                    DisplacementOperator.kernel,
                                    PolarisationOperator,
                                    parameters;
                                    linsolver = linsolver,
                                    nsteps = [nsteps, 1],
                                    FETypes = [FEType_D, FEType_P],
                                    target_residual = [tres, tres],
                                    solve_polarisation = polarisation,
                                    coupled = false,
                                    maxiterations = [maxits, 1])
        end
    else
        ## problem description
        PD = ProblemDescription("My problem")
        assign_unknown!(PD, u)
        if polarisation
            assign_unknown!(PD, V)
        end
        quadorder = (femorder-1)*2 + bonus_quadorder
        M = [Matrix(diagm(1.0 .+ ai)) for ai in eps0]
        kappa0 = 8.854e-12 # [[C^2/(N m^2)]]
        EO = EnergyOperator(M, MD.TensorC, polarisation ? MD.TensorE : nothing, [MD.data[j].kappar for j = 1:3] .* kappa0)
        if polarisation
            opid = assign_operator!(PD, NonlinearOperator(EO.eval_∂FW!, [grad(u), grad(V)]; sparse_jacobians = false, quadorder = quadorder, damping = damping))
            assign_operator!(PD, HomogeneousBoundaryData(V; regions = [4,5,6]))
        else
            opid = assign_operator!(PD, NonlinearOperator(EO.eval_∂FW!, [grad(u)]; sparse_jacobians = false, quadorder = quadorder, damping = damping))
        end
        assign_operator!(PD, HomogeneousBoundaryData(u; regions = [4,5,6]))

        FES = [FESpace{FEType_D}(xgrid), FESpace{FEType_P}(xgrid)]
        SC = nothing
        Solution = FEVector(FES; tags = PD.unknowns)
        last_residual = 0.0
        for j = 1 : nsteps
            println("Step $(j) out of $(nsteps)")
            M = [Matrix(diagm(1 .+ ai*j/nsteps)) for ai in eps0]
            update_M!(EO, M)
            if polarisation
                replace_operator!(PD, opid, NonlinearOperator(EO.eval_∂FW!, [grad(u), grad(V)]; sparse_jacobians = false, quadorder = quadorder))
            else
                replace_operator!(PD, opid, NonlinearOperator(EO.eval_∂FW!, [grad(u)]; sparse_jacobians = false, quadorder = quadorder))
            end
            Solution, SC = ExtendableFEM.solve(PD, FES, SC; init = Solution, maxiterations = maxits, return_config = true, method_linear = linsolver, target_residual = tres, damping = damping)
            last_residual = residual(SC)
        end
    end

    ## add computed data to dict
    resultd = deepcopy(d)
    resultd["eps0"] = eps0
    resultd["a"] = a
    resultd["solution"] = Solution
    resultd["residual"] = last_residual
    wsave(datadir(watson_datasubdir, filename), resultd)

    ###################
    ### POSTPROCESS ###
    ###################
    if d["postprocess"]
        d = postprocess(filename; Plotter = Plotter, cut_levels = geometry[4]/2, simple_cuts = true)
        resultd["angle"] = d["angle"]
        resultd["curvature"] = d["curvature"]
        wsave(datadir(watson_datasubdir, filename), resultd)
    else
        @info "skipping postprocessing... start it manually with: nanowire.postprocess(\"$filename\"; Plotter = PyPlot)"
    end

    return resultd
end


function get_lattice_misfit_nanowire(lcavg_case, MD::HeteroStructureData, geometry, full_nonlin)

    nregions = length(MD.data)
    lc_avg::Array{Float64,1} = zeros(Float64,3)
    eps0::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, nregions)
    r::Array{Float64,1} = geometry[1:3]
    a::Array{Float64,1} = zeros(Float64,3)
    E::Array{Float64,1} = zeros(Float64,3)

    lc = Array{Array{Float64,1}}(undef, nregions)
    for j = 1 : nregions
        lc[j] = MD.data[j].LatticeConstants
    end

    if nregions == 2
        # equalibrium strain for nanowire
        A_core = r[1]*r[2]
        A_stressor = r[1]*r[3]
        for j = 1 : 3
            #lc_avg[j] = (lc[1][j]*A_core + lc[2][j]*A_stressor)/(A_core + A_stressor)
            lc_avg[j] = 2*lc[1][j]*lc[2][j]/(lc[1][j] + lc[2][j])
        end
    elseif nregions == 3 # core and shell are a joint region
        # equilibrium strain for nanowire
        #A_core = 3*sqrt(3)/2 * r[1]^2
        #A_shell = 3*sqrt(3)/2 * (r[2]^2 + 2*r[1]*r[2])
        #A_stressor = sqrt(3)/2 * r[3] * (7*(r[1]+r[2]) + 3*r[3])
        d  = (r[1]+r[2])/sqrt(3)
        for region = 1 : nregions
            C = MD.TensorC[region].C
            E[region] = C[3,3] - 2 * C[1,3]^2/(C[1,1] + C[1,2])
        end

        for j = 1 : 3
            #lc_avg[j] = (lc[1][j]*A_core + lc[2][j]*A_shell + lc[3][j]*A_stressor)/(A_core + A_shell + A_stressor)
            lc_avg[j] = (lc[1][j]*(405*sqrt(3)*d^5*lc[3][j]*E[1]^3 - 54*d^2*r[3]*(d*(37*d + 12*sqrt(3)*r[3])*lc[1][j] -
                (81*d^2 + 30*sqrt(3)*d*r[3] + 16*r[3]^2)*lc[3][j])*E[1]^2*E[3] +
                6*d*r[3]^2*(-16*sqrt(3)*r[3]^2*(lc[1][j] - 6*lc[3][j]) + 72*d*r[3]*(-2*lc[1][j] + 5*lc[3][j]) +
                sqrt(3)*d^2*(-163*lc[1][j] + 314*lc[3][j]))*E[1]*E[3]^2 - 8*r[3]^3*(15*d^2 + 16*r[3]^2)*(lc[1][j] - 2*lc[3][j])*E[3]^3))/
                (lc[3][j]*(3*sqrt(3)*d*E[1] + 4*r[3]*E[3])*(135*d^4*E[1]^2 +
                12*d*r[3]*(17*sqrt(3)*d^2 + 27*d*r[3] + 8*sqrt(3)*r[3]^2)*E[1]*E[3] + 2*r[3]^2*(15*d^2 + 16*r[3]^2)*E[3]^2))
        end
    end

    for region = 1 : nregions
        eps0[region] = zeros(Float64,3)
        for j = 1 : 3
            if lcavg_case == 1
                a[j] = (lc[region][j] - lc[1][j])/lc[region][j]
            elseif lcavg_case == 2
                a[j] = (lc[region][j] - lc_avg[j])/lc[region][j]
            elseif lcavg_case == 3
                a[j] = (lc[region][j] - lc_avg[j])/lc_avg[j]
            end

            eps0[region][j] = a[j]
        end
    end

    return eps0, a
end

## load a result dictionary (to browse results in julia)
function load_data(d = nothing; watson_datasubdir = watson_datasubdir, kwargs...)
    ## complete parameter set
    if d === nothing
        d = get_defaults()
    end
    set_params!(d; kwargs...)

    # load dict from file
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    d = wload(datadir(watson_datasubdir, filename))
    return d
end

function export_vtk(d = nothing; upscaling = 0, kwargs...)
    if d === nothing
        d = load_data(d; kwargs...)
    end
    filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    @unpack solution, strainm, estrainm, eps0, polarisation = d
    repair_grid!(solution[1].FES.xgrid)
    if polarisation
        exportVTK(datadir(watson_datasubdir, filename_vtk), eps0, solution[1], solution[2]; EST = estrainm, strain_model = strainm, P0strain = true, upscaling = upscaling)
    else
        exportVTK(datadir(watson_datasubdir, filename_vtk), eps0, solution[1]; EST = estrainm, strain_model = strainm, P0strain = true, upscaling = upscaling)
    end
end

function postprocess(filename = nothing; watson_datasubdir = watson_datasubdir, Plotter = nothing, export_sol = true, cross_section_cuts = true, cut_levels = "auto", deform = true, simple_cuts = true, cut_npoints = 500, vol_cut = "auto", eps_gfind = 1e-12, upscaling = 0, kwargs...)

    if typeof(filename) <: Dict
        d = filename
        filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    elseif filename === nothing
        d = load_data(watson_datasubdir = watson_datasubdir; kwargs...)
        filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    else
        d = wload(datadir(watson_datasubdir, filename))
    end

    ## compute statistics
    @unpack solution, geometry, cross_section, rotate = d

    # cross_section_points = cross_section == 1 ? [0,-(geometry[1]+geometry[2])/sqrt(3)] : [0,-(geometry[1]+geometry[2])/2]
    # # rotate cross section points by angle theta
    # theta = rotate * pi/180
	# rotated_pointx = cos(theta)*cross_section_points[1] + sin(theta)*cross_section_points[2]
	# rotated_pointy = -sin(theta)*cross_section_points[1] + cos(theta)*cross_section_points[2]

    #bending_axis_end_points = [[rotated_pointx,rotated_pointy,0],[rotated_pointx,rotated_pointy,geometry[4]]]
    bending_axis_end_points = [[0,0,0],[0,0,geometry[4]]]
    angle, curvature, shifting_angle, ~ = compute_statistics(solution[1].FES.xgrid, solution[1], bending_axis_end_points, rotate)
    d["angle"] = angle
    d["curvature"] = curvature
    d["shifting_angle"] = shifting_angle

    # export vtk files
    if export_sol == true
        export_vtk(d)
    end

    ## save again
    @show filename
    wsave(datadir(watson_datasubdir, filename), d)

    ## compute cuts (only exported as vtu and png files)
    if cross_section_cuts
        @unpack strainm, nrefs = d
        if cut_levels == "auto"
            cut_levels = [0.3,0.4,0.5,0.6,0.7] * geometry[4]
        end
        if vol_cut == "auto"
            vol_cut = 4.0^(2-nrefs)
        end
        filename_cuts = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses) * "_CUTS/"
        mkpath(datadir(watson_datasubdir, filename_cuts))
        diam = geometry[1] + geometry[2]
        plane_points = [[-0.25*diam,-0.25*diam],[0.25*diam,-0.25*diam],[-0.25*diam,0.25*diam]] # = [X,Y] coordinates of the three points that define the cut plane
        # rescale cut_npoints based on geometry to fit cross section in a rectangle
        if cross_section == 1
            if rotate == 0 || rotate == 180
                cut_npoints = [cut_npoints, Int(ceil(cut_npoints*(2*diam/sqrt(3)+geometry[3])/diam))]
            elseif rotate == 90 || rotate == 270
                cut_npoints = [Int(ceil(cut_npoints*(2*diam/sqrt(3)+geometry[3])/diam)), cut_npoints]
            end
        elseif cross_section == 2
            if rotate == 0 || rotate == 180
                cut_npoints = [cut_npoints, Int(ceil(cut_npoints*sqrt(3)/2*(diam+geometry[3])/diam))]
            elseif rotate == 90 || rotate == 270
                cut_npoints = [Int(ceil(cut_npoints*sqrt(3)/2*(diam+geometry[3])/diam)), cut_npoints]
            end
        else
            cut_npoints = [cut_npoints, cut_npoints]
        end
        if simple_cuts # needs grid that triangulates cut_levels
            perform_simple_plane_cuts(datadir(watson_datasubdir, filename_cuts), solution, plane_points, cut_levels; eps0 = d["eps0"], eps_gfind = eps_gfind, deform = deform, cut_npoints = cut_npoints, only_localsearch = true, strain_model = strainm, Plotter = Plotter, upscaling = upscaling)
        else
            #perform_plane_cuts(datadir(watson_datasubdir, filename_cuts), solution, plane_points, cut_levels; strain_model = strainm, cut_npoints = cut_npoints, vol_cut = vol_cut, Plotter = Plotter)
        end
    end

    return d

end

end
