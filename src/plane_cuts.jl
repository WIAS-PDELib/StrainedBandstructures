using DrWatson

# modification from https://github.com/j-fu/GridVisualize.jl/blob/main/src/common.jl
function tet_x_plane!(ixcoord,ixvalues,iedges,pointlist,node_indices,planeq_values,function_values; tol=0.0)

    # If all nodes lie on one side of the plane, no intersection
    if (mapreduce(a->a< -tol,*,planeq_values) || mapreduce(a->a>tol,*,planeq_values))
        return 0
    end
    # Interpolate coordinates and function_values according to
    # evaluation of the plane equation
    edge_rule = local_celledgenodes(Tetrahedron3D)
    nxs=0
    n1::Int = 0
    n2::Int = 0
    for iedge = 1 : 6
        n1 = edge_rule[1,iedge]
        n2 = edge_rule[2,iedge]
        if planeq_values[n1]*planeq_values[n2]<tol
            nxs+=1
            t= planeq_values[n1]/(planeq_values[n1]-planeq_values[n2])
            for i=1:3
                ixcoord[i,nxs]=pointlist[i,node_indices[n1]]+t*(pointlist[i,node_indices[n2]]-pointlist[i,node_indices[n1]])
            end
            ixvalues[nxs]=function_values[node_indices[n1]]+t*(function_values[node_indices[n2]]-function_values[node_indices[n1]])
            # also remember the local edge numbers that had a succesful cut
            iedges[nxs] = iedge
        end
    end
    return nxs
end

# modification from marching_tetrahedra from https://github.com/j-fu/GridVisualize.jl/blob/main/src/common.jl
function plane_cut(xgrid::ExtendableGrid, plane_equation_coeffs; project_data = [], remesh::Bool = true, vol = 1, tol = 0.0)

    Tv=Float64
    Tp=Float64
    Tf=Int32

    all_ixfaces=Vector{Tf}(undef,0)
    all_ixbedges=Vector{Tf}(undef,0)
    all_ixcoord=Vector{Tp}(undef,0)
    all_ixvalues=Vector{Tv}(undef,0)

    planeq=zeros(4)
    ixcoord=zeros(3,6)
    ixvalues=zeros(6)
    ixbnd=zeros(Bool,6)
    iedges=zeros(Int,6)
    node_indices=zeros(Int32,4)

    plane_equation(plane,coord)=coord[1]*plane[1]+coord[2]*plane[2]+coord[3]*plane[3]+plane[4]

    xCoordinates = xgrid[Coordinates]
    xCellNodes=xgrid[CellNodes]
    nnodes = size(xCoordinates,2)
    func = zeros(Float64,1,nnodes)
    if project_data != []
        for FEB in project_data
            nodevalues!(func, FEB)
        end
    end

    # find nodes along boundary
    edge_on_boundary = zeros(Bool,num_sources(xgrid[GradientRobustMultiPhysics.EdgeNodes]))
    xFaceNodes = xgrid[FaceNodes]
    xFaceCells = xgrid[FaceCells]
    xFaceEdges = xgrid[FaceEdges]
    xCellEdges = xgrid[GradientRobustMultiPhysics.CellEdges]
    for j = 1 : size(xFaceNodes,2)
        if xFaceCells[2,j] == 0
            edge_on_boundary[xFaceEdges[:,j]] .= true
        end
    end

    function pushtris(ns,ixcoord,ixvalues,ixbnd)
        # number of intersection points can be 3 or 4
        if ns>=3
            last_i=length(all_ixvalues)
            for is=1:ns
                # todo: transform points onto z = 0 plane
                @views append!(all_ixcoord,ixcoord[:,is])
                push!(all_ixvalues,ixvalues[is])
            end
            if ixbnd[1] && ixbnd[2]
                append!(all_ixbedges,[last_i+1,last_i+2])
            end
            if ixbnd[2] && ixbnd[3]
                append!(all_ixbedges,[last_i+2,last_i+3])
            end
            if ixbnd[3] && ixbnd[1]
                append!(all_ixbedges,[last_i+3,last_i+1])
            end
            append!(all_ixfaces,[last_i+1,last_i+2,last_i+3])
            if ns==4
                append!(all_ixfaces,[last_i+3,last_i+2,last_i+4])
                if ixbnd[3] && ixbnd[2]
                    append!(all_ixbedges,[last_i+3,last_i+2])
                end
                if ixbnd[2] && ixbnd[4]
                    append!(all_ixbedges,[last_i+2,last_i+4])
                end
                if ixbnd[4] && ixbnd[3]
                    append!(all_ixbedges,[last_i+4,last_i+3])
                end
            end
        end
    end

    start_cell = 0
    for itet=1:size(xCellNodes,2)
        for i=1:4
            node_indices[i]=xCellNodes[i,itet]
        end
        
        @views map!(inode->plane_equation(plane_equation_coeffs,xCoordinates[:,inode]),planeq,node_indices)
        nxs=tet_x_plane!(ixcoord,ixvalues,iedges,xCoordinates,node_indices,planeq,func; tol=tol)
        if nxs >= 3
            start_cell = itet
        end
        # check if cutted edges are boundary edges
        for iedge = 1 : nxs
            ixbnd[iedge] = edge_on_boundary[xCellEdges[iedges[iedge],itet]]
        end
        # save triangles of cut
        pushtris(nxs,ixcoord,ixvalues,ixbnd)
    end

    ## reshape
    xCoordinates = reshape(all_ixcoord,3,Int(length(all_ixcoord)/3))
    nbedges = Int(length(all_ixbedges)/2)
    xBFaceNodes = reshape(all_ixbedges,2,nbedges)

    ## rotate normal around x axis such that n[2] = 0
    ## tan(alpha) = n[2]/n[3]
    alpha = atan(plane_equation_coeffs[2]/plane_equation_coeffs[3])
    #@info "rotating around x-axis with alpha = $alpha"
    plane_equation_coeffs[3] = sin(alpha)*plane_equation_coeffs[2] + cos(alpha)*plane_equation_coeffs[3]
    plane_equation_coeffs[2] = 0

    ## rotate normal around y axis such that n[1] = 0
    ## tan(alpha) = n[2]/n[3]
    beta = -atan(plane_equation_coeffs[1]/plane_equation_coeffs[3])
   # @info "rotating around x-axis with beta = $beta"
    plane_equation_coeffs[3] = -sin(beta)*plane_equation_coeffs[1] + cos(beta)*plane_equation_coeffs[3]
    plane_equation_coeffs[1] = 0

    ## rotate coordinates
    oldcoords = zeros(Float64,3)
    R =  [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)] * [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]
    for j = 1 : size(xCoordinates,2)
        for k = 1 : 3
            oldcoords[k] = xCoordinates[k,j]
            xCoordinates[k,j] = 0
        end
        for n = 1 : 3, k = 1 : 3
            xCoordinates[n,j] += R[n,k] * oldcoords[k]
        end
    end

    ## restrict coordinates to [x,y] plane
    z = xCoordinates[end,3] # are all the same
    xCoordinates = xCoordinates[1:2,:]

    if remesh
        # problem: boundary nodes are not unique, there are nodes with the same coordinates that we have
        # the following lines look for unique boundary nodes and remaps the BFaceNodes
        bnodes = unique(xBFaceNodes[:])
        node_remap = zeros(Int,length(xCoordinates))
        xCoordinatesB = zeros(Float64,0)
        already_known::Int = 0
        nbnodes::Int = 0
        for node in bnodes
            already_known = 0
            for knode = 1 : nbnodes
                if abs((xCoordinatesB[2*knode-1] - xCoordinates[1,node]).^2 + (xCoordinatesB[2*knode] - xCoordinates[2,node]).^2) < 1e-16
                    already_known = knode
                    break
                end
            end
            if already_known == 0
                append!(xCoordinatesB,xCoordinates[:,node])
                nbnodes += 1
                node_remap[node] = nbnodes
            else
                node_remap[node] = already_known
            end
        end
        @views xBFaceNodes[:] .= node_remap[xBFaceNodes[:]]
        xCoordinatesB = reshape(xCoordinatesB,2,Int(length(xCoordinatesB)/2))
        CenterPoint = sum(xCoordinatesB, dims = 2) ./ size(xCoordinatesB,2)

        # call Grid generator
        xgrid = SimplexGridFactory.simplexgrid(Triangulate;
            points=xCoordinatesB,
            bfaces=xBFaceNodes,
            bfaceregions=ones(Int32,nbedges),
            regionpoints=CenterPoint',
            regionnumbers=[1],
            regionvolumes=[vol])
    else
        ncells = Int(length(all_ixfaces)/3)
        xCellNodes = reshape(all_ixfaces,3,ncells)
        xgrid=ExtendableGrid{Float64,Int32}()
        xgrid[Coordinates] = xCoordinates
        xgrid[CellNodes] =  xCellNodes
        xgrid[CellGeometries] = VectorOfConstants(Triangle2D,ncells);
        xgrid[CellRegions]=ones(Int32,ncells)
        xgrid[BFaceRegions]=ones(Int32,nbedges)
        xgrid[BFaceNodes]=xBFaceNodes
        xgrid[BFaceGeometries]=VectorOfConstants(Edge1D,nbedges)
        xgrid[CoordinateSystem]=Cartesian2D
    end
    
    return xgrid, R, z, start_cell
end



## computes two grids: one boundary conforming Delaunay grid for finite volume methods
## and one uniform non-boundary-conforming one that comprises the cut;
## also a function is returned that transforms 3D coordinates to 2D coordinates on the cut
function get_cutgrids(xgrid, plane_equation_coeffs; npoints = 100, vol_cut = 1.0)
    @time cut_grid, R, z, start_cell = plane_cut(xgrid, plane_equation_coeffs; remesh = true, vol = vol_cut)

    ## define regular grid on cut plane where polarisation should be interpolated
    xmin = minimum(view(cut_grid[Coordinates],1,:))
    xmax = maximum(view(cut_grid[Coordinates],1,:))
    ymin = minimum(view(cut_grid[Coordinates],2,:))
    ymax = maximum(view(cut_grid[Coordinates],2,:))
    @info "Creating uniform grid for bounding box ($xmin,$xmax) x ($ymin,$ymax)"
    h_uni = [xmax - xmin,ymax-ymin] ./ (npoints-1)
    xgrid_uni = simplexgrid(collect(xmin:h_uni[1]:xmax),collect(ymin:h_uni[2]:ymax))

    # mapping from 3D coordinates on cut_grid to 2D coordinates on xgrid_uni
    invR::Matrix{Float64} = inv(R)
    function xtrafo!(x3D,x2D)
        for j = 1 : 3
            x3D[j] = invR[j,1] * x2D[1] + invR[j,2] * x2D[2] + invR[j,3] * z  # invR * [x2D[1],x2D[2],z]
        end
        return nothing
    end

    return cut_grid, xgrid_uni, xtrafo!, start_cell
end


function perform_simple_plane_cuts(target_folder_cut, Solution_original, plane_points, cut_levels; 
    strain_model = NonlinearStrain3D, 
    export_uniform_data = true,
    cut_npoints = [200,200],
    only_localsearch = true,
    eps_gfind = 1e-11,
    deform = false,
    eps0 = nothing,
    EST = AnisotropicDiagonalPrestrain,
    cut_direction = 3, # 1 = y-z-plane, 2 = x-z-plane, 3 = x-y-plane (default)
    Plotter = nothing,
    upscaling = 0,
    tight_box = true,
    do_simplecut_plots = Plotter !== nothing,
    do_uniformcut_plots = Plotter !== nothing)

    if cut_direction == 1
        direction_string = "x"
    elseif cut_direction == 2
        direction_string = "y"
    elseif cut_direction == 3
        direction_string = "z"
    else 
        @error "cut_direction needs to be 1,2 or 3"
    end

    PolarisationGradient = nothing
    if upscaling > 0
        xgrid_original = Solution_original[1].FES.xgrid
        xgrid = uniform_refine(xgrid_original, upscaling)
        FESup = FESpace{eltype(Solution_original[1].FES)}(xgrid)
        Solution = FEVector(FESup)
        @info "INTERPOLATING TO UPSCALED GRID..."
        interpolate!(Solution[1], Solution_original[1])
        #interpolate!(Solution[2], Solution_original[2])
        DisplacementGradient = continuify(Solution[1], Gradient)
        @info "STARTING CUTTING..."
        if length(Solution) > 1
            PolarisationGradient = continuify(Solution[2], Gradient)
        end
    else
        Solution = Solution_original
        DisplacementGradient = continuify(Solution[1], Gradient)
        if length(Solution) > 1
            PolarisationGradient = continuify(Solution[2], Gradient)
        end
        xgrid = Solution_original[1].FES.xgrid
    end


    gridplot(xgrid; Plotter=Plotter)
    
    ### first find faces that lie in cut_levels

    # faces with all nodes with the same z-coordinate are assigned to a z-level (otherwise z-level -1 is set)
    xCoordinates = xgrid[Coordinates]
    xCellRegions = xgrid[CellRegions]
    xFaceNodes = xgrid[FaceNodes]
    xFaceCells = xgrid[FaceCells]
    nfaces::Int = num_sources(xFaceNodes)
    z4Faces = zeros(Float64,nfaces)
    nnodes4face = max_num_targets_per_source(xFaceNodes)
    nnodes = size(xCoordinates,2)
    component_names = ["XX","YY","ZZ","YZ","XZ","XY"]

    z::Float64 = 0
    for face = 1 : nfaces
        z = xCoordinates[cut_direction,xFaceNodes[1,face]]
        for k = 2 : nnodes4face
            if abs(xCoordinates[cut_direction,xFaceNodes[k,face]] - z) > eps_gfind
                z = -1
                break;
            end
        end
        z4Faces[face] = z
    end

    ## these are the levels that have faces
    unique_levels = sort(unique(z4Faces))[2:end]
    @info "These are the possible $(direction_string) levels for simple cuts:"
    for l = 1 : length(unique_levels)
        if unique_levels[l] in cut_levels
            print(" * ")
        else
            print("   ")
        end
        println("level = $(unique_levels[l]) | nfaces on level = $(length(findall(abs.(z4Faces .- unique_levels[l]) .< eps_gfind)))")
    end

    ## find three points on the plane z = cut_level and evaluate displacement at points of plane
    xref = [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    cells = zeros(Int,3)
    PE = PointEvaluator([id(1)], Solution)
    CF = CellFinder(xgrid)
    plane_equation_coeffs = zeros(Float64,4)

    strain = zeros(Float64,6)
    gradient = zeros(Float64,9)


    ## interpolate eps0 on original grid
    eps0_fefunc_orig = FEVector(FESpace{L2P0{3}}(xgrid))
    ncells_orig = num_cells(xgrid)
    for cell = 1 : ncells_orig
        if typeof(eps0[xCellRegions[cell]]) <: Real
            for k = 1 : 3
                eps0_fefunc_orig.entries[(cell-1)*3 + k] = eps0[xCellRegions[cell]]
            end
        else
            eps0_fefunc_orig.entries[(cell-1)*3 + 1] = eps0[xCellRegions[cell]][1]
            eps0_fefunc_orig.entries[(cell-1)*3 + 2] = eps0[xCellRegions[cell]][2]
            eps0_fefunc_orig.entries[(cell-1)*3 + 3] = eps0[xCellRegions[cell]][3]
        end
    end

    function get_cutgrid(nodes4level, faces4level)
        ## map original 3D coordinates onto 2D plane
        nnodes_cut = length(nodes4level)
        xCoordinatesCut3D = zeros(Float64,3,nnodes_cut)
        for j = 1 : nnodes_cut
            for k = 1 : 3
                xCoordinatesCut3D[k,j] = xCoordinates[k,nodes4level[j]]
            end
        end

        node_permute = zeros(Int32,size(xCoordinates,2))
        node_permute[nodes4level] = 1 : length(nodes4level)

        ## construct simple cut grid
        cut_grid=ExtendableGrid{Float64,Int32}()
        cut_grid[Coordinates] = xCoordinatesCut3D
        cut_grid[CellNodes] = node_permute[xFaceNodes[:,faces4level]] # view not possible here
        cut_grid[CellGeometries] = VectorOfConstants{ElementGeometries,Int}(Triangle2D,length(faces4level));
        cut_grid[CellRegions] = xgrid[CellRegions][xgrid[FaceCells][1,faces4level]]
        cut_grid[CoordinateSystem] = Cartesian2D
        #cut_grid[BFaceRegions] = ones(Int32,0)
        #cut_grid[BFaceNodes] = zeros(Int32,2,0)
        #cut_grid[BFaceCells] = zeros(Int32,2,0)
        cut_grid[BFaceRegions] = ones(Int32,1)
        cut_grid[BFaceNodes] = Matrix{Int32}([1 2;])
        #cut_grid[BFaceCells] = zeros(Int32,2,0)
        cut_grid[BFaceGeometries] = VectorOfConstants{ElementGeometries,Int}(Edge1D, 0)
        return cut_grid
    end

    for l = 1 : length(cut_levels)

        cut_level = cut_levels[l]
        @info "ENTERING cut_level = $(cut_level)"

        # make cut_level subdirectory
        mkpath(datadir(target_folder_cut, "z=$(cut_level)/"))
        target_folder_cut_level = target_folder_cut * "z=$(cut_level)/"

        ## find faces for this cut_level
        faces4level = findall(abs.(z4Faces .- cut_level) .< eps_gfind)
        if length(faces4level) < 1
            @warn "found no faces in grid that lie in level = $(cut_level), skipping this cut level..."
            continue
        end
        @show length(faces4level)
        nodes4level = unique(view(xFaceNodes,:,faces4level))
        nnodes_cut = length(nodes4level)
        start_cell::Int = xFaceCells[1,faces4level[1]]

        ## define plane equation coefficients and rotation to map 3D coordinates on cut plane to 2D coordinates

        # 1:3 = normal vector
        #   4 = - normal vector ⋅ point on plane
        # find normal vector of displaced plane defined by the three points x[1], x[2] and x[3] 
        if cut_direction == 1
            x = [[cut_level, plane_points[1][1],plane_points[1][2]],[cut_level,plane_points[2][1],plane_points[2][2]],[cut_level,plane_points[3][1],plane_points[3][2]]]
            a = 2
            b = 3
            c = 1
        elseif cut_direction == 2
            x = [[plane_points[1][1],cut_level,plane_points[1][2]],[plane_points[2][1],cut_level,plane_points[2][2]],[plane_points[3][1],cut_level,plane_points[3][2]]]
            a = 3
            b = 1
            c = 2
        elseif cut_direction == 3
            x = [[plane_points[1][1],plane_points[1][2],cut_level],[plane_points[2][1],plane_points[2][2],cut_level],[plane_points[3][1],plane_points[3][2],cut_level]]
            a = 1
            b = 2
            c = 3
        end
        result = deepcopy(x[1])

        if deform
            ## the three points that define the plane on the reference state
            ## are displaced to get three points on the transformed plane
            for i = 1 : 3
                # find cell
                cells[i] = gFindLocal!(xref[i], CF, x[i]; icellstart = start_cell, eps = eps_gfind)
                @info cells[i]
                if cells[i] == 0
                    @warn "local search cell search unexpectedly failed, using brute force..."
                    cells[i] = gFindBruteForce!(xref[i], CF, x[i])
                end
                @assert cells[i] > 0
                # evaluate displacement
                evaluate_bary!(result,PE,xref[i],cells[i])
                ## displace point
                x[i] .+= result
            end
        end

        ## normal = (x[1]-x[2]) × (x[1]-x[3])
        plane_equation_coeffs[1]  = (x[1][2] - x[2][2]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[1] -= (x[1][3] - x[2][3]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[2]  = (x[1][3] - x[2][3]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[2] -= (x[1][1] - x[2][1]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[3]  = (x[1][1] - x[2][1]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[3] -= (x[1][2] - x[2][2]) * (x[1][1] - x[3][1])
        plane_equation_coeffs ./= sqrt(sum(plane_equation_coeffs.^2))
        plane_equation_coeffs[4] = -sum(x[1] .* plane_equation_coeffs[1:3])

        ## rotate normal around axis a such that n[b] = 0
        ## tan(alpha) = n[b]/n[c]
        alpha = atan(plane_equation_coeffs[b]/plane_equation_coeffs[c])
        @info "rotating around axis $a with alpha = $alpha"
        plane_equation_coeffs[c] = sin(alpha)*plane_equation_coeffs[b] + cos(alpha)*plane_equation_coeffs[c]
        plane_equation_coeffs[b] = 0

        ## rotate normal around axis b such that n[a] = 0
        ## tan(alpha) = n[a]/n[c]
        beta = -atan(plane_equation_coeffs[a]/plane_equation_coeffs[c])
        @info "rotating around axis $b with beta = $beta"
        plane_equation_coeffs[c] = -sin(beta)*plane_equation_coeffs[a] + cos(beta)*plane_equation_coeffs[c]
        plane_equation_coeffs[a] = 0

        if cut_direction == 1
            R = [cos(beta) 0 sin(beta); -sin(beta) 0 cos(beta); 0 1 0] * [cos(alpha) -sin(alpha) 0; 0 0 1; sin(alpha) cos(alpha) 0]
        elseif cut_direction == 2
            R = [ 0 0 1; cos(beta) sin(beta) 0; -sin(beta) cos(beta) 0] * [0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha); 1 0 0]
        elseif cut_direction == 3
            R = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)] * [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]
        end

        ## that is the rotation matrix that maps points from the reference domain to the x-y plane times a fixed z coordinates
        @show R


        #### SIMPLE CUT GRIDS ###

        ## interpolate identity at nodes4level
        nodevals = nodevalues(Solution[1], Identity; continuous = true, nodes = nodes4level)
        if length(Solution) > 1
            nodevals_P = nodevalues(Solution[2], Identity; continuous = true, nodes = nodes4level)
        else
            nodevals_P = nothing
        end

        ## get displaced cut_grid from function defined above
        cut_grid = get_cutgrid(nodes4level, faces4level)

        ## get subgrid for each region
        subgrid1 = subgrid(cut_grid, [1,2]; boundary = false)
        subgrid2 = subgrid(cut_grid, [3]; boundary = false)

        ## get parent nodes for each subgrid
        subnodes1 = subgrid1[ExtendableGrids.NodeParents]
        subnodes2 = subgrid2[ExtendableGrids.NodeParents]
        subgrid1[Coordinates] = cut_grid[Coordinates][:,subnodes1]
        subgrid2[Coordinates] = cut_grid[Coordinates][:,subnodes2]

        ## interpolate data on cut_grid
        @info "Interpolating data on cut mesh..."
        nodevals_gradient = nodevalues(Solution[1], Gradient; nodes = nodes4level)
        if length(Solution) > 1
            nodevals_E = nodevalues(Solution[2], Gradient; nodes = nodes4level)
        else
            nodevals_E = nothing
        end
        nodevals_ϵu = zeros(Float64,6,nnodes_cut)
        nodevals_ϵu_elastic = zeros(Float64,6,nnodes_cut)
        
        ## compute nodevalues for nodes of each subgrid
        nodevals_gradient1 = nodevalues(Solution[1], Gradient; regions = [1,2], nodes = nodes4level[subnodes1])
        nodevals_gradient2 = nodevalues(Solution[1], Gradient; regions = [3], nodes = nodes4level[subnodes2])
        nodevals_ϵu1 = zeros(Float64,6,size(nodevals_gradient1,2))
        nodevals_ϵu2 = zeros(Float64,6,size(nodevals_gradient2,2))

        ## compute nodevalues for polarisatin/electri field on subgrids
        if length(Solution) > 1
            nodevals_P1 = nodevalues(Solution[2], Identity; regions = [1,2], nodes = nodes4level[subnodes1])
            nodevals_P2 = nodevalues(Solution[2], Identity; regions = [3], nodes = nodes4level[subnodes2])
            nodevals_E1 = nodevalues(Solution[2], Gradient; regions = [1,2], nodes = nodes4level[subnodes1])
            nodevals_E2 = nodevalues(Solution[2], Gradient; regions = [3], nodes = nodes4level[subnodes2])
        end

        ## now displace the grid if deform is true
        xCoordinatesCut3D = cut_grid[Coordinates]
        if deform
            for j = 1 : nnodes_cut
                for k = 1 : 3
                    xCoordinatesCut3D[k,j] += nodevals[k,j]
                end
            end
            subgrid1[Coordinates] = cut_grid[Coordinates][:,subnodes1]
            subgrid2[Coordinates] = cut_grid[Coordinates][:,subnodes2]
        end

        ## interpolate eps0 on cut_grid
        eps0_fefunc = FEVector(FESpace{L2P0{3}}(cut_grid))
        ncells_cut = num_cells(cut_grid)
        xCellRegions_cut = cut_grid[CellRegions]
        for cell = 1 : ncells_cut
            if typeof(eps0[xCellRegions_cut[cell]]) <: Real
                for k = 1 : 3
                    eps0_fefunc.entries[(cell-1)*3 + k] = eps0[xCellRegions_cut[cell]]
                end
            else
                eps0_fefunc.entries[(cell-1)*3 + 1] = eps0[xCellRegions_cut[cell]][1]
                eps0_fefunc.entries[(cell-1)*3 + 2] = eps0[xCellRegions_cut[cell]][2]
                eps0_fefunc.entries[(cell-1)*3 + 3] = eps0[xCellRegions_cut[cell]][3]
            end
        end

        ## calculate strain from gradient interpolation on (undisplaced) cut
        for nv in [[nodevals_gradient,nodevals_ϵu,0], [nodevals_gradient1,nodevals_ϵu1,1], [nodevals_gradient2,nodevals_ϵu2,3]]
            nv_∇u = nv[1]
            nv_ϵu = nv[2]
            region = nv[3]
            for j = 1 : size(nv_∇u,2)
                eval_strain!(strain, view(nv_∇u,:,j), strain_model)
                if region !== 0
                    eval_elastic_strain!(strain, eps0[region], EST)
                end
                nv_ϵu[:,j] .= strain
                for k = 4 : 6
                    nv_ϵu[k,j] /= 2
                end
            end
        end

        ## calculate jumps along interface of core and stressor
        nodes_at_interface = intersect(subnodes1, subnodes2)

        ## choose coordinate that changes sign (quick and dirty criterion)
        xinterface = cut_grid[Coordinates][1, nodes_at_interface]
        label = "jumps along x"
        if prod(extrema(xinterface)) > 0
            xinterface = cut_grid[Coordinates][2, nodes_at_interface]
        label = "jumps along y"
        end

        ## compute jumps
        jumps = zeros(Float64, 6, length(nodes_at_interface))
        for j = 1 : length(nodes_at_interface)
            n = nodes_at_interface[j]
            j1 = findfirst(==(n), subnodes1)
            j2 = findfirst(==(n), subnodes2)
            jumps[:,j] =  nodevals_ϵu1[:,j1] .- nodevals_ϵu2[:,j2]
        end

        ## sort jumps w.r.t to coordinates
        P = sortperm(xinterface)
        xinterface = xinterface[P]
        jumps = jumps[:,P]

        ## report
        for c = 1 : 6
            @info "minimal/maximal jump at interface in $(component_names[c]) component: $(extrema(jumps[c,:]))"
        end
        
        if Plotter !== nothing
            interfacegrid = simplexgrid(xinterface)
            jump_fefunc = FEVector(FESpace{H1P1{1}}(interfacegrid), entries = jumps[1,:])

            plt = GridVisualizer(; Plotter = Plotter, layout = (1,1), clear = true, size = (800,600))
            colors = [:black, :blue, :red, :green, :yellow, :magenta]
            plt1 = nothing
            for c = 1 : 6
                jump_fefunc.entries .= jumps[c,:]
                plt1 = scalarplot!(plt[1,1], interfacegrid, jump_fefunc[1]; clear = false, label = component_names[c], color = colors[c], markershape = :circle, markevery = 1, markersize = 10, legend = :cb)
            end
            filename = target_folder_cut_level * "jumps_along_interface.png"
            if isdefined(Plotter,:savefig)
                Plotter.savefig(filename)
            else
                GridVisualize.save(filename, plt1; Plotter = Plotter)
            end
        end

        ## get 2D coordinates of the simple grid by applying the rotation R
        cut_grid2D = deepcopy(cut_grid)
        xCoordinatesCutPlane = zeros(Float64,2,nnodes_cut)
        zs = zeros(Float64, nnodes_cut)
        for j = 1 : nnodes_cut
            for k = 1 : 3
                xCoordinatesCutPlane[1,j] += R[a,k] * (xCoordinatesCut3D[k,j])
                xCoordinatesCutPlane[2,j] += R[b,k] * (xCoordinatesCut3D[k,j])
                zs[j] += R[cut_direction,k] * (xCoordinatesCut3D[k,j])
            end
        end
        z = sum(zs) / nnodes_cut
        cut_grid2D[Coordinates] = xCoordinatesCutPlane
        CF2D = CellFinder(cut_grid2D)

        ## write data into csv file
        @info "Exporting cut data for cut_level = $(cut_level)..."
        kwargs = Dict()
        kwargs1 = Dict()
        kwargs2 = Dict()
        kwargs[:cellregions] = cut_grid[CellRegions]
        kwargs[:displacement] = nodevals
        kwargs[:grad_displacement] = nodevals_gradient
        kwargs[:strain] = nodevals_ϵu
 #       kwargs[:elastic_strain] = nodevals_ϵu_elastic
        if length(Solution) > 1
            kwargs[:polarisation] = nodevals_P
            kwargs[:electric_field] = nodevals_E
            kwargs1[:electric_field] = nodevals_E1
            kwargs2[:electric_field] = nodevals_E2
        end

        kwargs1[:elastic_strain] = nodevals_ϵu1
        kwargs2[:elastic_strain] = nodevals_ϵu2
        ExtendableGrids.writeVTK(target_folder_cut_level * "simple_cut_$(cut_level)_data" * (deform ? "_deformed.vtu" : ".vtu"), cut_grid; kwargs...)
        ExtendableGrids.writeVTK(target_folder_cut_level * "simple_cut_$(cut_level)_data_subgrid1" * (deform ? "_deformed.vtu" : ".vtu"), subgrid1; kwargs1...)
        ExtendableGrids.writeVTK(target_folder_cut_level * "simple_cut_$(cut_level)_data_subgrid2" * (deform ? "_deformed.vtu" : ".vtu"), subgrid2; kwargs2...)


        ## data into txt files
        subgrid1[Coordinates] = xCoordinatesCutPlane[:,subnodes1]
        subgrid2[Coordinates] = xCoordinatesCutPlane[:,subnodes2]
        coords1 = subgrid1[Coordinates]
        coords2 = subgrid2[Coordinates]
        @info "Writing coordinates of subgrids..."
        filename_eAB = target_folder_cut_level * "coordinates_subgrid1.dat"
        io = open(filename_eAB, "w")
        for n = 1 : length(subnodes1)
            @printf(io, "%.6f %.6f\n",coords1[1,n],coords1[2,n])
        end
        close(io)
        filename_eAB = target_folder_cut_level * "coordinates_subgrid2.dat"
        io = open(filename_eAB, "w")
        for n = 1 : length(subnodes2)
            @printf(io, "%.6f %.6f\n",coords2[1,n],coords2[2,n])
        end
        close(io)
        for c = 1 : 6
            @info "Writing elastic strain distribution file for e_elastic$(component_names[c]) on subgrids..."
            filename_eAB = target_folder_cut_level * 'e' * component_names[c] * "_elastic_subgrid1.dat"
            io = open(filename_eAB, "w")
            for n = 1 : length(subnodes1)
                @printf(io, "%.6f\n",nodevals_ϵu1[c,n])
            end
            close(io)
            filename_eAB = target_folder_cut_level * 'e' * component_names[c] * "_elastic_subgrid2.dat"
            io = open(filename_eAB, "w")
            for n = 1 : length(subnodes2)
                @printf(io, "%.6f\n",nodevals_ϵu2[c,n])
            end
            close(io)
        end

        if length(Solution) > 1
            ## write polarisation potential and eletric field maps into txt files
            @info "Writing polarisation potential file for e$(component_names[c]) on subgrids..."
            filename_P = target_folder_cut_level * "P_subgrid1.dat"
            io = open(filename_P, "w")
            for n = 1 : length(subnodes1)
                @printf(io, "%.6f\n",nodevals_P1[n])
            end
            close(io)
            filename_P = target_folder_cut_level * "P_subgrid2.dat"
            io = open(filename_P, "w")
            for n = 1 : length(subnodes2)
                @printf(io, "%.6f\n",nodevals_P2[n])
            end
            close(io)

            pcomponents = ["Ex","Ey","Ez"]
            for c = 1 : 3
                @info "Writing electric field file for $(pcomponents[c]) on subgrids..."
                filename_Ec = target_folder_cut_level * pcomponents[c] * "_subgrid1.dat"
                io = open(filename_Ec, "w")
                #@printf(io, "%s\n", component_names[c])
                for n = 1 : length(subnodes1)
                    @printf(io, "%.6f\n",nodevals_E1[c,n])
                end
                close(io)
                filename_Ec = target_folder_cut_level * pcomponents[c] * "_subgrid2.dat"
                io = open(filename_Ec, "w")
                #@printf(io, "%s\n", component_names[c])
                for n = 1 : length(subnodes2)
                    @printf(io, "%.6f\n",nodevals_E2[c,n])
                end
                close(io)
            end
        end

        xmin = minimum(view(xCoordinatesCutPlane,a,:))
        xmax = maximum(view(xCoordinatesCutPlane,a,:))
        ymin = minimum(view(xCoordinatesCutPlane,b,:))
        ymax = maximum(view(xCoordinatesCutPlane,b,:))
        if do_simplecut_plots
            @info "Plotting data on simple cut grid..."
            labels = ["ux","uy","uz"]
            for j = 1 : 3
                plt = scalarplot(cut_grid2D, view(nodevals,j,:), Plotter = Plotter; xlimits = (xmin-2,xmax+2), ylimits = (ymin-2,ymax+2), title = "$(labels[j]) on cut", fignumber = 1)
                filename = target_folder_cut_level * "simple_cut_$(cut_level)_$(labels[j]).png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
	                GridVisualize.save(filename, plt; Plotter = Plotter)
                end
            end
            for k = 1 : 6
                plt = scalarplot(cut_grid2D, view(nodevals_ϵu,k,:), Plotter = Plotter; xlimits = (xmin-2,xmax+2), ylimits = (ymin-2,ymax+2), title = "ϵ_$(component_names[k]) on cut", fignumber = 1)
                filename = target_folder_cut_level * "simple_cut_$(cut_level)_ϵ$(component_names[k]).png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
	                GridVisualize.save(filename, plt; Plotter = Plotter)
                end
                plt = scalarplot([subgrid1,subgrid2], cut_grid2D, [view(nodevals_ϵu1,k,:),view(nodevals_ϵu2,k,:)], Plotter = Plotter; title = "ϵ_elastic_$(component_names[k]) on cut", fignumber = 1)
                filename = target_folder_cut_level * "simple_cut_subgrid_$(cut_level)_ϵ_elastic_$(component_names[k]).png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
	                GridVisualize.save(filename, plt; Plotter = Plotter)
                end
                # scalarplot(cut_grid2D, view(nodevals_ϵu_elastic,k,:), Plotter = Plotter; title = "ϵ_elastic_$(component_names[k]) on cut", fignumber = 1)
                # if isdefined(Plotter,:savefig)
                #     Plotter.savefig(target_folder_cut_level * "simple_cut_$(cut_level)_ϵ_elastic_$(component_names[k]).png")
                # end
            end
        end


        #### UNIFORM CUT GRIDS ###
        if export_uniform_data

            if deform == false
                ## start with the same grid coordinates from the simple grid
                ## (possibly deformed if deform = true)
                xCoordinatesCutPlane = xCoordinatesCut3D
            else
                xCoordinatesCutPlane = xCoordinatesCutPlane
            end

            ## restrict coordinates to [x,y] plane
            xmin = minimum(view(xCoordinatesCutPlane,a,:))
            xmax = maximum(view(xCoordinatesCutPlane,a,:))
            ymin = minimum(view(xCoordinatesCutPlane,b,:))
            ymax = maximum(view(xCoordinatesCutPlane,b,:))

            ## define bounding box and uniform cut grid
            ## In case tight_box == false then include some marginal cells at the boundary
            d = [xmax - xmin,ymax - ymin]
            if tight_box == false
                xmin -= d[1]*0.01
                xmax += d[1]*0.01
                ymin -= d[2]*0.01
                ymax += d[2]*0.01
            end
            hx_uni = d ./ (cut_npoints[1]-1)
            hy_uni = d ./ (cut_npoints[2]-1)
            Xuni = zeros(Float64,0)
            Yuni = zeros(Float64,0)
            for j = 0 : cut_npoints[1] - 1
                push!(Xuni, xmin + hx_uni[1]*j)
            end
            for j = 0 : cut_npoints[2] - 1
                push!(Yuni, ymin + hy_uni[2]*j)
            end
            @info "Creating uniform grid with $hx_uni x $hy_uni points for bounding box ($xmin,$xmax) x ($ymin,$ymax)"
            xgrid_uni = simplexgrid(Xuni,Yuni)
            xCoordinatesUni = xgrid_uni[Coordinates]
            nnodes_uni = size(xCoordinatesUni,2)
            if deform == false
                xgrid_uni[Coordinates] = [xCoordinatesUni; z*ones(Float64, nnodes_uni)']
            else
                ## how to get the z coordinates for the cut plane (also of the nodes that are outside)?
            end

            ## interpolate data on uniform cut_grid
            @info "Interpolating data on uniform cut mesh..."
            FES2D = FESpace{H1P1{3}}(xgrid_uni)
            FES2D_∇u = FESpace{H1P1{9}}(xgrid_uni)
            FES2D_ϵ = FESpace{H1P1{6}}(xgrid_uni)
            FES2D_P = FESpace{H1P1{1}}(xgrid_uni)
            FES2D_∇P = FESpace{H1P1{3}}(xgrid_uni)
            CutSolution_u = FEVector(FES2D)
            CutSolution_∇u = FEVector(FES2D_∇u)
            CutSolution_ϵu = FEVector(FES2D_ϵ)
            CutSolution_ϵu_elastic = FEVector(FES2D_ϵ)
            CutSolution_P = FEVector(FES2D_P)
            CutSolution_∇P = FEVector(FES2D_∇P)
            eps0_fefunc_uni = FEVector(FESpace{H1P1{3}}(xgrid_uni))

            # mapping from 3D coordinates on cut_grid to 2D coordinates on xgrid_uni
            invR::Matrix{Float64} = inv(R)
            z_offset = 0.0
            function xtrafo!(x3D,x2D)
                for j = 1 : 3
                    x3D[j] = invR[j,1] * x2D[1] + invR[j,2] * x2D[2] + invR[j,3] * z
                end
                return nothing
            end

            # mapping from 3D coordinates on simple cut_grid to 2D coordinates on xgrid_uni
            function xtrafo2!(x2D_out,x2D_in)
                x2D_out[1] = x2D_in[1] + xmin
                x2D_out[2] = x2D_in[2] + ymin
                return nothing
            end

            #interpolate ϵ0 onto uniform grid
            if deform
                displace_mesh!(xgrid, Solution[1]; magnify = 1)
                lazy_interpolate!(eps0_fefunc_uni[1], eps0_fefunc_orig, [id(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                lazy_interpolate!(CutSolution_u[1], Solution, [id(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                lazy_interpolate!(CutSolution_∇u[1], DisplacementGradient, [id(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                if length(Solution) > 1
                    lazy_interpolate!(CutSolution_P[1], Solution, [id(2)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                    lazy_interpolate!(CutSolution_∇P[1], PolarisationGradient, [id(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                end
                displace_mesh!(xgrid, Solution[1]; magnify = -1)
            else
                lazy_interpolate!(eps0_fefunc_uni[1], eps0_fefunc_orig, [id(1)]; not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                lazy_interpolate!(CutSolution_u[1], Solution, [id(1)]; not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                lazy_interpolate!(CutSolution_∇u[1], DisplacementGradient, [id(1)]; not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                if length(Solution) > 1
                    lazy_interpolate!(CutSolution_P[1], Solution, [id(2)]; not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                    lazy_interpolate!(CutSolution_∇P[1], PolarisationGradient, [id(1)]; not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
                end
            end

            ## postprocess gradient to gradients on undisplaced mesh
            ## and calculate strain values
            for j = 1 : nnodes_uni
                for k = 1 : 9
                    gradient[k] = CutSolution_∇u.entries[(k-1)*nnodes_uni + j]
                end

                ## phi = x + u(x) 
                ## dphi/dx = I + grad_x(u) =: F
                ## grad_phi(u) we have evaluated above
                ## 
                ## want: grad_x(u(x)) = grad_phi(u(phi)) * grad_x(phi) = grad_phi(u) * (I + grad_x(u))

                ## ==> grad_x(u) (I - grad_phi(u)) = grad_phi(u)
                ## ==> grad_x(u) = grad_phi(u) * inv(I - grad_phi(u))

                ## inv(F)*grad_x(u) = inv(F)*(F - I) = I - inv(F) = grad_phi(u)

                ## inv(I - grad_phi(u)) - I = grad_x(u)

                ###
                #gradmatrix = (-inv([gradient[1]-1 gradient[2] gradient[3];
                #              gradient[4] gradient[5]-1 gradient[6];
                #              gradient[7] gradient[8] gradient[9]-1]) - [1 0 0;0 1 0; 0 0 1])'
                #for k = 1 : 9
                #    CutSolution_∇u.entries[(k-1)*nnodes_uni + j] = gradmatrix[k]
                #end
                eval_strain!(strain, gradient, strain_model)
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes_uni + j] = strain[k]
                end
                eval_elastic_strain!(strain, eps0_fefunc_uni.entries[[j, nnodes_uni + j, 2*nnodes_uni + j]], EST)
                for k = 1 : 6
                    CutSolution_ϵu_elastic.entries[(k-1)*nnodes_uni + j] = strain[k]
                end
                for k = 4 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes_uni + j] /= 2
                    CutSolution_ϵu_elastic.entries[(k-1)*nnodes_uni + j] /= 2
                end
            end

            ## write material map into txt files
            @info "Writing material map file..."
            filename_material = target_folder_cut_level * "map.dat"
            io = open(filename_material, "w")
            xCellNodesUni = xgrid_uni[CellNodes]
            nnodes_uni = size(xCoordinatesUni,2)
            xdim = size(xCoordinatesUni,1)
            @printf(io, "%d %d\n", cut_npoints[1], cut_npoints[2])
            @printf(io, "CORE STRESSOR OUTSIDE\n")
            cell::Int = 0
            region::Int = 0
            x3D = zeros(Float64,3)
            xCellRegionsUniform = zeros(Int32,num_cells(xgrid_uni))
            xCellRegionsSimpleCut = cut_grid[CellRegions]

            ## get 2D coordinates of the uniform grid by applying the inverse of rotation R
            if deform == false
                xCoordinatesCut3D = xgrid_uni[Coordinates]
                xCoordinatesCut2D = zeros(Float64,2,nnodes_uni)
                invR = inv(R[[a,b,c],:])
                for j = 1 : nnodes_uni
                    z = 0
                    for k = 1 : 3
                        xCoordinatesCut2D[1,j] += invR[1,k] * (xCoordinatesCut3D[k,j])
                        xCoordinatesCut2D[2,j] += invR[2,k] * (xCoordinatesCut3D[k,j])
                        z += invR[3,k] * (xCoordinatesCut3D[k,j])
                    end
                end
            else
                xCoordinatesCut2D = xCoordinatesUni
            end

            xtest = zeros(Float64,2)
            for cu = 1 : num_cells(xgrid_uni)
                ## compute center
                fill!(xtest,0)
                for j = 1 : 2, n = 1 : 3
                    xtest[j] += xCoordinatesCut2D[j, xCellNodesUni[n,cu]]
                end
                xtest ./= 3.0
                #xtest[1] += xmin
                #xtest[2] += ymin
                cell = gFindLocal!(xref[1], CF2D, xtest; trybrute = true, eps = eps_gfind)
                xCellRegionsUniform[cu] = cell == 0 ? 0 : xCellRegionsSimpleCut[cell]
            end
            xgrid_uni[CellRegions] = xCellRegionsUniform

            for n = 1 : nnodes_uni
                for j = 1 : 2
                    xtest[j] = xCoordinatesCut2D[j,n]
                end
                #xtest[1] += xmin
                #xtest[2] += ymin
                cell = gFindLocal!(xref[1], CF2D, xtest; trybrute = true, eps = eps_gfind)
                region = cell == 0 ? 0 : xCellRegionsSimpleCut[cell]
                if region == 1
                    @printf(io, "1 0 0 ")
                elseif region == 2
                    @printf(io, "1 0 0 ")
                elseif region == 3
                    @printf(io, "0 1 0 ")
                elseif region == 0
                    @printf(io, "0 0 1 ")
                end
                #for j = 1 : xdim
                #    @printf(io, "%.6f ",xCoordinatesUni[j,n])
                #end
                @printf(io, "\n")
            end
            close(io)

            ## write data into vtk file
            @info "Writing data into vtk file..."
            kwargs = Dict()
            kwargs[:cellregions] = xgrid_uni[CellRegions]
            kwargs[:displacement] = view(nodevalues(CutSolution_u[1], Identity),:,:)
            kwargs[:grad_displacement] = view(nodevalues(CutSolution_∇u[1], Identity),:,:)
            kwargs[:strain] = view(nodevalues(CutSolution_ϵu[1], Identity),:,:)
            kwargs[:elastic_strain] = view(nodevalues(CutSolution_ϵu_elastic[1], Identity),:,:)
            if length(Solution) > 1
                kwargs[:polarisation] = view(nodevalues(CutSolution_P[1], Identity),:,:)
                kwargs[:electric_field] = view(nodevalues(CutSolution_∇P[1], Identity),:,:)
            end
            ExtendableGrids.writeVTK(target_folder_cut_level * "uniform_cut_$(cut_level)_data" * (deform ? "_deformed.vtu" : ".vtu"), xgrid_uni; kwargs...)
           

            ## replacing NaN with 1e30 so that min/max calculation works
            replace!(CutSolution_u.entries, NaN=>1e30)
            replace!(CutSolution_∇u.entries, NaN=>1e30)
            replace!(CutSolution_ϵu.entries, NaN=>1e30)
            replace!(CutSolution_ϵu_elastic.entries, NaN=>1e30)
            replace!(CutSolution_P.entries, NaN=>1e30)
            replace!(CutSolution_∇P.entries, NaN=>1e30)

            ## write strain distribution maps into txt files
            for c = 1 : 6
                @info "Writing strain distribution file for e$(component_names[c])..."
                filename_eAB = target_folder_cut_level * 'e' * component_names[c] * ".dat"
                io = open(filename_eAB, "w")
                #@printf(io, "%s\n", component_names[c])
                for n = 1 : nnodes_uni
                    @printf(io, "%.6f\n",CutSolution_ϵu.entries[(c-1)*nnodes_uni+n])
                end
                close(io)
                @info "Writing elastic strain distribution file for e_elastic$(component_names[c])..."
                filename_eAB = target_folder_cut_level * 'e' * component_names[c] * "_elastic.dat"
                io = open(filename_eAB, "w")
                #@printf(io, "%s\n", component_names[c])
                for n = 1 : nnodes_uni
                    @printf(io, "%.6f\n",CutSolution_ϵu_elastic.entries[(c-1)*nnodes_uni+n])
                end
                close(io)
            end

            ## write polarisation potential and eletric field maps into txt files
            @info "Writing polarisation potential file for e$(component_names[c])..."
            filename_P = target_folder_cut_level * "P.dat"
            io = open(filename_P, "w")
            for n = 1 : nnodes_uni
                @printf(io, "%.6f\n",CutSolution_P.entries[n])
            end
            close(io)
            pcomponents = ["Ex","Ey","Ez"]
            for c = 1 : 3
                @info "Writing electric field file for $(pcomponents[c])..."
                filename_Ec = target_folder_cut_level * pcomponents[c] * ".dat"
                io = open(filename_Ec, "w")
                #@printf(io, "%s\n", component_names[c])
                for n = 1 : nnodes_uni
                    @printf(io, "%.6f\n",CutSolution_∇P.entries[(c-1)*nnodes_uni+n])
                end
                close(io)
            end

            ## plot displacement, strain and polarisation on uniform cut grid
            if do_uniformcut_plots

                xgrid_uni2D = deepcopy(xgrid_uni)
                xgrid_uni2D[Coordinates] = xCoordinatesCut2D
                #gridplot(xgrid_uni2D; Plotter = Plotter)

                @info "Plotting data on uniform cut grid..."
                uxmin::Float64 = 1e30
                uxmax::Float64 = -1e30
                uymin::Float64 = 1e30
                uymax::Float64 = -1e30
                uzmin::Float64 = 1e30
                uzmax::Float64 = -1e30
                Pmin::Float64 = 1e30
                Pmax::Float64 = -1e30
                ϵmax = -1e30*ones(Float64,6)
                ϵmin = 1e30*ones(Float64,6)
                ϵmax_elastic = -1e30*ones(Float64,6)
                ϵmin_elastic = 1e30*ones(Float64,6)
                nnodes_uni = size(xgrid_uni[Coordinates],2)
                for j = 1 : nnodes_uni
                    if abs(CutSolution_u.entries[j]) < 1e10
                        if length(Solution) > 1
                            Pmin = min(Pmin,CutSolution_P[1][j])
                            Pmax = max(Pmax,CutSolution_P[1][j])
                        end
                        uxmin = min(uxmin,CutSolution_u[1][j])
                        uymin = min(uymin,CutSolution_u[1][nnodes_uni+j])
                        uzmin = min(uzmin,CutSolution_u[1][2*nnodes_uni+j])
                        uxmax = max(uxmax,CutSolution_u[1][j])
                        uymax = max(uymax,CutSolution_u[1][nnodes_uni+j])
                        uzmax = max(uzmax,CutSolution_u[1][2*nnodes_uni+j])
                        if abs(CutSolution_ϵu.entries[j]) < 1e10
                            for k = 1 : 6
                                ϵmax[k] = max(ϵmax[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                                ϵmin[k] = min(ϵmin[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                                ϵmax_elastic[k] = max(ϵmax_elastic[k],CutSolution_ϵu_elastic[1][(k-1)*nnodes_uni+j])
                                ϵmin_elastic[k] = min(ϵmin_elastic[k],CutSolution_ϵu_elastic[1][(k-1)*nnodes_uni+j])
                            end
                        end
                    end
                end
                @show uxmin, uxmax
                plt = scalarplot(xgrid_uni2D, view(CutSolution_u.entries,1:nnodes_uni), Plotter = Plotter; flimits = (uxmin,uxmax), title = "ux on cut", fignumber = 1)
                filename = target_folder_cut_level * "uniform_cut_$(cut_level)_ux.png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
                    GridVisualize.save(filename, plt; Plotter = Plotter)
                end
                plt = scalarplot(xgrid_uni2D, view(CutSolution_u.entries,nnodes_uni+1:2*nnodes_uni), Plotter = Plotter; flimits = (uymin,uymax), title = "uy on cut", fignumber = 1)
                filename = target_folder_cut_level * "uniform_cut_$(cut_level)_uy.png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
                    GridVisualize.save(filename, plt; Plotter = Plotter)
                end
                plt = scalarplot(xgrid_uni2D, view(CutSolution_u.entries,2*nnodes_uni+1:3*nnodes_uni), Plotter = Plotter; flimits = (uzmin,uzmax), title = "uz on cut", fignumber = 1)
                filename = target_folder_cut_level * "uniform_cut_$(cut_level)_uz.png"
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(filename)
                else
                    GridVisualize.save(filename, plt; Plotter = Plotter)
                end
                if length(Solution) > 1
                    plt = scalarplot(xgrid_uni2D, CutSolution_P.entries, Plotter = Plotter; flimits = (Pmin,Pmax), title = "Polarisation on cut", fignumber = 1)
                    filename = target_folder_cut_level * "uniform_cut_$(cut_level)_P.png"
                    if isdefined(Plotter,:savefig)
                        Plotter.savefig(filename)
                    else
                        GridVisualize.save(filename, plt; Plotter = Plotter)
                    end
                end
                for k = 1 : 6
                    plt = scalarplot(xgrid_uni2D, view(CutSolution_ϵu.entries,(k-1)*nnodes_uni+1:k*nnodes_uni), Plotter = Plotter; flimits = (ϵmin[k],ϵmax[k]), title = "ϵ_$(component_names[k]) on cut", fignumber = 1)
                    filename = target_folder_cut_level * "uniform_cut_$(cut_level)_ϵ$(component_names[k]).png"
                    if isdefined(Plotter,:savefig)
                        Plotter.savefig(filename)
                    else
                        GridVisualize.save(filename, plt; Plotter = Plotter)
                    end
                end
                for k = 1 : 6
                    plt = scalarplot(xgrid_uni2D, view(CutSolution_ϵu_elastic.entries,(k-1)*nnodes_uni+1:k*nnodes_uni), Plotter = Plotter; flimits = (ϵmin_elastic[k],ϵmax_elastic[k]), title = "ϵ_$(component_names[k]) on cut", fignumber = 1)
                    filename = target_folder_cut_level * "uniform_cut_$(cut_level)_ϵ_elastic$(component_names[k]).png"
                    if isdefined(Plotter,:savefig)
                        Plotter.savefig(filename)
                    else
                        GridVisualize.save(filename, plt; Plotter = Plotter)
                    end
                end
            end
        end
    end
end


function perform_plane_cuts(target_folder_cut, Solution, plane_points, cut_levels; strain_model = NonlinearStrain3D, cut_npoints = 100, 
    only_localsearch = true, vol_cut = 16, eps_gfind = 1e-11, Plotter = nothing)

    xgrid = Solution[1].FES.xgrid

    ## find three points on the plane z = cut_level and evaluate displacement at points of plane
    @info "Calculating coefficients of plane equations for cuts at levels $(cut_levels)"
    xref = [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    cells = zeros(Int,3)
    PE = PointEvaluator(Solution[1], Identity)
    CF = CellFinder(xgrid)

    plane_equation_coeffs = Array{Array{Float64,1},1}(undef, length(cut_levels))
    for l = 1 : length(cut_levels)

        ## define plane equation coefficients
        ## 1:3 = normal vector
        ##   4 = - normal vector ⋅ point on plane
        ## find normal vector of displaced plane defined by the three points x[1], x[2] and x[3] 
        cut_level = cut_levels[l]

        x = [[plane_points[1][1],plane_points[1][2],cut_level],[plane_points[2][1],plane_points[2][2],cut_level],[plane_points[3][1],plane_points[3][2],cut_level]]

        result = deepcopy(x[1])
        for i = 1 : 3
            # find cell
            cells[i] = gFindLocal!(xref[i], CF, x[i]; icellstart = 1, eps = eps_gfind)
            if cells[i] == 0
                cells[i] = gFindBruteForce!(xref[i], CF, x[i])
            end
            @assert cells[i] > 0
            # evaluate displacement
            evaluate!(result,PE,xref[i],cells[i])
            ## displace point
            x[i] .+= result
        end

        plane_equation_coeffs[l] = zeros(Float64,4)
        plane_equation_coeffs[l][1]  = (x[1][2] - x[2][2]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[l][1] -= (x[1][3] - x[2][3]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[l][2]  = (x[1][3] - x[2][3]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[l][2] -= (x[1][1] - x[2][1]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[l][3]  = (x[1][1] - x[2][1]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[l][3] -= (x[1][2] - x[2][2]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[l] ./= sqrt(sum(plane_equation_coeffs[l].^2))
        plane_equation_coeffs[l][4] = -sum(x[1] .* plane_equation_coeffs[l][1:3])
    end

    ## displace grid
    displace_mesh!(xgrid, Solution[1])

    ## cut displaced grid at plane
    for l = 1 : length(cut_levels)
        cut_level = cut_levels[l]
        @info "Cutting domain at z = $cut_level with plane equation coefficients $(plane_equation_coeffs[l])"
        @time cut_grid, xgrid_uni, xtrafo!, start_cell = get_cutgrids(xgrid, plane_equation_coeffs[l]; npoints = cut_npoints, vol_cut = vol_cut)

        # plot boundary-conforming Delaunay cut mesh (suitable for FV)
       # @info "Plotting Delaunay cut mesh..."
        # gridplot(cut_grid, Plotter = Plotter, title = "Delaunay mesh of cut", fignumber = 1)

        ## interpolate data on uniform cut_grid
        @info "Interpolating data on uniform cut mesh..."
        FES2D = FESpace{H1P1{3}}(xgrid_uni)
        FES2D_∇u = FESpace{H1P1{9}}(xgrid_uni)
        FES2D_ϵ = FESpace{H1P1{6}}(xgrid_uni)
        FES2D_P = FESpace{H1P1{1}}(xgrid_uni)
        CutSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D)
        CutSolution_∇u = FEVector{Float64}("∇u (on 2D cut at z = $(cut_level))", FES2D_∇u)
        CutSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_ϵ)
        CutSolution_P = FEVector{Float64}("P (on 2D cut at z = $(cut_level))", FES2D_P)
        @time lazy_interpolate!(CutSolution_u[1], Solution, [id(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
        @time lazy_interpolate!(CutSolution_∇u[1], Solution, [grad(1)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
          
        ## calculate strain from gradient interpolation on cut
        nodevals = nodevalues(CutSolution_∇u[1], Identity)
        strain = zeros(Float64,6)
        nnodes = size(nodevals,2)
        for j = 1 : nnodes
            if nodevals[1,j] == NaN
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes+j] = NaN
                end
            else
                eval_strain!(strain,view(nodevals,:,j), strain_model)
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes+j] = strain[k]
                end
            end
        end
        if length(Solution) > 1
            @time lazy_interpolate!(CutSolution_P[1], Solution, [id(2)]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
        end

        ## write data into csv file
        @info "Writing data into csv file..."
        writeVTK!(target_folder_cut * "cut_$(cut_level)_data.vtu", [CutSolution_u[1],CutSolution_∇u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity, Identity])
        writeCSV!(target_folder_cut * "cut_$(cut_level)_data.txt", [CutSolution_u[1],CutSolution_∇u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity, Identity], seperator = "\t")

        ## replacing NaN with 1e30 so that min/max calculation works
        replace!(CutSolution_u.entries, NaN=>1e30)
        replace!(CutSolution_∇u.entries, NaN=>1e30)
        replace!(CutSolution_ϵu.entries, NaN=>1e30)
        replace!(CutSolution_P.entries, NaN=>1e30)

        ## plot displacement, strain and polarisation on uniform cut grid
        @info "Plotting data on uniform cut grid..."
        uxmin::Float64 = 1e30
        uxmax::Float64 = -1e30
        uymin::Float64 = 1e30
        uymax::Float64 = -1e30
        uzmin::Float64 = 1e30
        uzmax::Float64 = -1e30
        Pmin::Float64 = 1e30
        Pmax::Float64 = -1e30
        ϵmax = -1e30*ones(Float64,6)
        ϵmin = 1e30*ones(Float64,6)
        nnodes_uni = size(xgrid_uni[Coordinates],2)
        for j = 1 : nnodes_uni
            if abs(CutSolution_u.entries[j]) < 1e10
                if length(Solution) > 1
                    Pmin = min(Pmin,CutSolution_P[1][j])
                    Pmax = max(Pmax,CutSolution_P[1][j])
                end
                uxmin = min(uxmin,CutSolution_u[1][j])
                uymin = min(uymin,CutSolution_u[1][nnodes_uni+j])
                uzmin = min(uzmin,CutSolution_u[1][2*nnodes_uni+j])
                uxmax = max(uxmax,CutSolution_u[1][j])
                uymax = max(uymax,CutSolution_u[1][nnodes_uni+j])
                uzmax = max(uzmax,CutSolution_u[1][2*nnodes_uni+j])
                if abs(CutSolution_ϵu.entries[j]) < 1e10
                    for k = 1 : 6
                        ϵmax[k] = max(ϵmax[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                        ϵmin[k] = min(ϵmin[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                    end
                end
            end
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,1:nnodes_uni), Plotter = Plotter; flimits = (uxmin,uxmax), title = "ux on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_ux.png")
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,nnodes_uni+1:2*nnodes_uni), Plotter = Plotter; flimits = (uymin,uymax), title = "uy on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_uy.png")
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,2*nnodes_uni+1:3*nnodes_uni), Plotter = Plotter; flimits = (uzmin,uzmax), title = "uz on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_uz.png")
        end
        if length(Solution) > 1
            scalarplot(xgrid_uni, CutSolution_P.entries, Plotter = Plotter; flimits = (Pmin,Pmax), title = "Polarisation on cut", fignumber = 1)
            if isdefined(Plotter,:savefig)
                Plotter.savefig(target_folder_cut * "cut_$(cut_level)_P.png")
            end
        end
        for k = 1 : 6
            scalarplot(xgrid_uni, view(CutSolution_ϵu.entries,(k-1)*nnodes_uni+1:k*nnodes_uni), Plotter = Plotter; flimits = (ϵmin[k],ϵmax[k]), title = "ϵu[$k] on cut", fignumber = 1)
            if isdefined(Plotter,:savefig)
                Plotter.savefig(target_folder_cut * "cut_$(cut_level)_ϵ$k.png")
            end
        end
    end
end
