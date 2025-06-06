function solve_lowlevel(
    xgrid::ExtendableGrid,                          # grid
    boundary_operator,                                   # Operator for boundary data (Dirichlet)
    displacement_operator,                          # Operator for displacement
    polarisation_operator,                          # Operator for polarisation
    emb_params;                                     # embedding parameters (operators must depend on them!)
    periodic_boundary_operator = nothing,             # operator for periodic boundary
    FETypes = [H1P1{size(xgrid[Coordinates],1)}],   # FETypes (default: P1)
    linsolver = ExtendableSparse.MKLPardiso,
    nsteps = ones(Int,length(FETypes)),             # number of embedding steps (parameters are scaled by nsteps equidistant steps within 0:1)
    damping = 0,                                    # damping in Newton iteration
    target_residual = 1e-12*ones(Int,length(FETypes)),
    solve_polarisation = false,                     # also solve for polarisation ?
    coupled = false,                                # solve coupled system in each nonlinear iteration ? (not needed as long as polarisation does not influence displacement)
    bonus_quadorder = 2,
    use_sparsity_detection = false,                 # sparsity detection for jacobians
    maxiterations = 20*ones(Int,length(FETypes)),
    fixed_penalty = 1e60
    )

    ## generate FES spaces
    Tv = Float64
    Ti = Int32
    FES = Array{FESpace{Tv,Ti},1}(undef, 0)
    push!(FES, FESpace{FETypes[1]}(xgrid))

    solve_polarisation = (length(FETypes) == 2) * solve_polarisation
    coupled = coupled * solve_polarisation

    if solve_polarisation
        push!(FES, FESpace{FETypes[2]}(xgrid))
    end

    ## generate FEVectors and FEMatrix
    SolutionOld = FEVector(FES)     # storage for previous solution
    Solution::FEVector{Float64} = FEVector(FES)        # storage for current and final solution
    if coupled
        S = FEMatrix(FES)       # system matrix
        rhs = FEVector(FES)     # system rhs
    else
        S = FEMatrix(FES[1])       # system matrix
        rhs = FEVector(FES[1])     # system rhs
    end

    @info "Starting low level solver...
            ndofs = $(length(Solution.entries))
            polarisation = $solve_polarisation (coupled = $coupled)"

    ## prepare linear solver
	LP = LinearProblem(S.entries.cscmatrix, rhs.entries)
    linsolve = nothing

    ## generate quadrature rule
    EG = xgrid[UniqueCellGeometries][1]
    order = get_polynomialorder(FETypes[1], EG)
    qf = QuadratureRule{Float64, EG}(2*(order-1) + bonus_quadorder)
    weights::Vector{Float64} = qf.w
    nweights::Int = length(weights)
    cellvolumes::Array{Float64,1} = xgrid[CellVolumes]
    cellregions::Array{Int,1} = xgrid[CellRegions]

    ## prepare dofmaps and FE evaluators
    ndofs_u::Int = get_ndofs(ON_CELLS, FETypes[1], EG)      # displacement
    ∇u::FEEvaluator{Float64} = FEEvaluator(FES[1], Gradient, qf)
    vals_∇u::Array{Float64,3} = ∇u.cvals
    celldofs_u::Matrix{Int32} = Matrix(FES[1][CellDofs])
    indices_∇u = 1:size(vals_∇u, 1)

    if solve_polarisation # polarisation
        ndofs_P::Int = get_ndofs(ON_CELLS, FETypes[2], EG)  
        ∇P::FEEvaluator{Float64} = FEEvaluator(FES[2], Gradient, qf)
        vals_∇P::Array{Float64,3} = ∇P.cvals
        celldofs_P::Matrix{Int32} = Matrix(FES[2][CellDofs])
        indices_∇P = 1:size(vals_∇P, 1)
    end

    ## local matrix and vector structures
    AUU::Array{Float64,2} = zeros(Float64, ndofs_u, ndofs_u)
    bU::Array{Float64,1} = zeros(Float64, ndofs_u)
    SE::ExtendableSparseMatrix{Float64,Int64} = S.entries
    SolutionE::Array{Float64,1} = Solution.entries

    offsets::Array{Int,1} = [0, FES[1].ndofs]
    if solve_polarisation
        push!(offsets, FES[1].ndofs + FES[2].ndofs)
    end

    ## prepare automatic differentation of displacement operator
    dim = size(xgrid[Coordinates],1)
    if solve_polarisation
        @assert dim == polarisation_operator.dim
    end
    argsizes_displacement = [dim^2,dim^2,Int(dim^2+dim*(dim+1)/2)]
    argsizes_polarisation = [dim,dim+dim^2,Int(dim+dim^2+dim*(dim+1)/2)]
    QP = QPInfos(xgrid)
    displacement_operator4qp = (result,input) -> displacement_operator(result,input,QP)
    polarisation_operator4qp = (result,input) -> polarisation_operator(result,input,QP)
    result_displacement::Array{Float64,1} = Vector{Float64}(undef,argsizes_displacement[1])
    input_displacement::Array{Float64,1} = Vector{Float64}(undef,argsizes_displacement[3])
    Dresult_displacement = DiffResults.JacobianResult(result_displacement,input_displacement)
    result_polarisation::Array{Float64,1} = Vector{Float64}(undef,argsizes_polarisation[1])
    input_polarisation::Array{Float64,1} = Vector{Float64}(undef,argsizes_polarisation[3])
    Dresult_polarisation = DiffResults.JacobianResult(result_polarisation,input_polarisation)

    if (use_sparsity_detection)
        sparsity_pattern = Symbolics.jacobian_sparsity(displacement_operator4qp,result_displacement,input_displacement)
        jac = Float64.(sparse(sparsity_pattern))
        colors = matrix_colors(jac) 
        cfg_displacement = ForwardColorJacCache(displacement_operator4qp,input_displacement,nothing;
                dx = nothing,
                colorvec = colors,
                sparsity = nothing)
    else
        cfg_displacement = ForwardDiff.JacobianConfig(displacement_operator4qp, result_displacement, input_displacement, ForwardDiff.Chunk{argsizes_displacement[3]}())
        cfg_polarisation = ForwardDiff.JacobianConfig(polarisation_operator4qp, result_polarisation, input_polarisation, ForwardDiff.Chunk{argsizes_polarisation[3]}())
    end
    jac_displacement::Array{Float64,2} = DiffResults.jacobian(Dresult_displacement)
    value_displacement::Array{Float64,1} = DiffResults.value(Dresult_displacement)
    jac_polarisation::Array{Float64,2} = DiffResults.jacobian(Dresult_polarisation)
    value_polarisation::Array{Float64,1} = DiffResults.value(Dresult_polarisation)
    
    ## ASSEMBLY LOOP
    bdofs = nothing
    bdofsP = nothing
    eval_∇u::Array{Float64,1} = zeros(Float64, argsizes_displacement[3])
    eval_∇P::Array{Float64,1} = zeros(Float64, argsizes_polarisation[3])
    jac_displacement_view = view(jac_displacement,indices_∇u,indices_∇u)
    ncells::Int = num_cells(xgrid)
    tempV::Array{Float64,1} = zeros(Float64, argsizes_displacement[1])
    tempP::Array{Float64,1} = zeros(Float64, argsizes_polarisation[1])

    function update_displacement(; coupled = false)

        fill!(S[1,1],0)
        fill!(rhs[1],0)
        if coupled
            fill!(S[1,2],0)
        end

        for cell = 1 : ncells

            ## update FE basis evaluators
            QP.item = cell
            QP.region = cellregions[cell]
            ∇u.citem[] = cell
            update_basis!(∇u) 
            dofs_u = view(celldofs_u, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u in current solution at current quadrature point
                fill!(eval_∇u, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇u[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_displacement, displacement_operator4qp, eval_∇u, cfg_displacement)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_displacement,displacement_operator4qp,result_displacement,eval_∇u,cfg_displacement)
                end

                # update matrix
                for j = 1 : ndofs_u
                    # multiply with jacobian
                    mul!(tempV, jac_displacement_view, view(vals_∇u,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_u
                        AUU[j,k] += dot(tempV, view(vals_∇u,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempV, jac_displacement_view, view(eval_∇u, indices_∇u))
                displacement_operator(value_displacement, eval_∇u, QP)
                tempV .-= value_displacement
                for j = 1 : ndofs_u
                    bU[j] += dot(tempV, view(vals_∇u,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_u
                for k = 1 : ndofs_u
                    if abs(AUU[j,k]) > 0
                        rawupdateindex!(SE, +, AUU[j,k] * cellvolumes[cell] , dofs_u[k], dofs_u[j])
                    end
                end
                rhs.entries[dofs_u[j]] += bU[j] * cellvolumes[cell]
            end

            fill!(AUU, 0)
            fill!(bU, 0)
        end

        flush!(SE)

        ## apply changes for periodic boundary
        if length(periodic_boundary_operator) > 0
            for op in periodic_boundary_operator
                ExtendableFEM.build_assembler!(op, [Solution[1], Solution[1]])
                op.assembler(SE, rhs.entries, true, true)
            end
        end

        ## apply penalties for boundary condition
        if boundary_operator !== nothing
            assemble!(boundary_operator, FES[1])
            bdofs = fixed_dofs(boundary_operator)
            for dof in bdofs
                SE[dof,dof] = fixed_penalty
                rhs.entries[dof] = 0
                Solution.entries[dof] = 0
            end
            flush!(SE)
        end
        
    end

    if solve_polarisation
        jac_polarisation_view = view(jac_polarisation,indices_∇P,dim^2+1:dim^2+dim)
        APP::Array{Float64,2} = zeros(Float64, ndofs_P, ndofs_P)
        bP::Array{Float64,1} = zeros(Float64, ndofs_P)
    end

    function update_polarisation(; coupled = false)

        if coupled
            offset = FES[1].ndofs
            fill!(S[2,1],0)
            fill!(S[2,2],0)
            fill!(rhs[2],0)
        else
            offset = 0
            fill!(S[1,1],0)
            fill!(rhs[1],0)
        end

        for cell = 1 : ncells

            ## update FE basis evaluators
            item_information[1] = cell
            item_information[2] = cell
            item_information[3] = cellregions[cell]
            pola_op = polarisation_operator4region(item_information)
            ∇u.citem[] = cell
            update_basis!(∇u) 
            ∇P.citem[] = cell
            update_basis!(∇P) 
            dofs_u = view(celldofs_u, : , cell)
            dofs_P = view(celldofs_P, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u and ∇P in current solution at current quadrature point
                fill!(eval_∇P, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇P[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end
                for j = 1 : ndofs_P, k in indices_∇P
                    eval_∇P[9+k] += SolutionE[dofs_P[j]+offsets[2]] * vals_∇P[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_polarisation, pola_op, eval_∇P, cfg_polarisation)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_polarisation,pola_op,result_polarisation,eval_∇P,cfg_polarisation)
                end

                # update matrix
                for j = 1 : ndofs_P
                    # multiply with jacobian
                    mul!(tempP, jac_polarisation_view, view(vals_∇P,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_P
                        APP[j,k] += dot(tempP, view(vals_∇P,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempP, jac_polarisation_view, view(eval_∇P, indices_∇P))
                pola_op(value_polarisation, eval_∇P)
                tempP .-= value_polarisation
                for j = 1 : ndofs_P
                    bP[j] += dot(tempP, view(vals_∇P,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_P
                for k = 1 : ndofs_P
                    if abs(APP[j,k]) > 0
                        rawupdateindex!(SE, +, APP[j,k] * cellvolumes[cell] , dofs_P[k], dofs_P[j])
                    end
                end
                rhs.entries[dofs_P[j]] += bP[j] * cellvolumes[cell]
            end

            fill!(APP, 0)
            fill!(bP, 0)
        end

        ## apply penalties for boundary condition
        fixed_dofsP = [offset + 1]
        for dof in fixed_dofsP
            SE[dof,dof] = fixed_penalty
            rhs.entries[dof] = 0 # fixed_penalty * Solution.entries[dof]
            Solution.entries[dof] = 0
        end

        flush!(SE)
    end

    ## prepare other stuff for loop
    emb_params_target = deepcopy(emb_params)
    factorization = nothing
    residual = zeros(Float64, size(SE,1))
    nlres::Float64 = 0
    linres::Float64 = 0
    change::Float64 = 0
    time_assembly::Float64 = 0
    time_solver::Float64 = 0
    iteration::Int = 0

    for step = 1 : nsteps[1]
        emb_params .= nsteps[1] == 1 ? emb_params_target : (step-1)/(nsteps[1]-1) .* emb_params_target

        @info "Entering embedding step = $step/$(nsteps[1]), emb_params = $emb_params"
            
        iteration = 0
        time_solver = 0
        while iteration <= maxiterations[1]
            # update system
            time_assembly = @elapsed update_displacement(; coupled = coupled)
            if iteration == 0 && step == 1
                linsolve = init(LP, linsolver)
            end
        
            # compute nonlinear residual
            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            if boundary_operator !== nothing
                view(residual, bdofs) .= 0
            end
            if length(periodic_boundary_operator) > 0
                for op in periodic_boundary_operator
                    view(residual, fixed_dofs(op)) .= 0
                end
            end

            nlres = norm(residual)
            @info "         ---> nonlinear iteration $iteration, linres = $linres, nlres = $nlres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
            if nlres < target_residual[1]
                @info "target reached, stopping..."
                break
            end

            # update factorization and solve
            if damping > 0
                SolutionOld.entries .= Solution.entries
            end

            time_solver = @elapsed begin
				linsolve.A = S.entries.cscmatrix
				linsolve.b = rhs.entries
				x = LinearSolve.solve!(linsolve)
                Solution.entries .= x.u
            end
            iteration += 1

            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            if boundary_operator !== nothing
                view(residual, bdofs) .= 0
            end
            if length(periodic_boundary_operator) > 0
                for op in periodic_boundary_operator
                    view(residual, fixed_dofs(op)) .= 0
                end
            end
            linres = norm(residual)

            if damping > 0
                view(Solution[1]) .= (1 - damping) * Solution.entries + damping * SolutionOld.entries
                change = norm(SolutionOld.entries - Solution.entries)
            end
        end
    end

    if !coupled && solve_polarisation
        @info "Entering decoupled polarisation solve..."
        S = FEMatrix(FES[2])       # system matrix
        rhs = FEVector(FES[2])     # system rhs
        SE = S.entries
        residual = zeros(Float64, size(SE,1))

        time_assembly = @elapsed update_polarisation(; coupled = coupled)
        factorization=linsolver(SE)
        factorize!(factorization, SE)

        time_solver = @elapsed begin
            ExtendableSparse.update!(factorization)
            ldiv!(view(Solution[2]), factorization, rhs.entries)
        end
        iteration += 1

        fill!(residual, 0)
        mul!(residual, SE, view(Solution[2]))
        residual .-= rhs.entries
        view(residual, fixed_dofsP) .= 0
        linres = norm(residual)

        @info "         ---> linres = $linres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
        
    end

    return Solution, nlres

end



function solve_lowlevel_parallel(
    xgrid::ExtendableGrid,                          # grid
    boundary_operator,                                   # Operator for boundary data (Dirichlet)
    displacement_operator,                          # Operator for displacement
    polarisation_operator,                          # Operator for polarisation
    emb_params;                                     # embedding parameters (operators must depend on them!)
    periodic_boundary_operator = nothing,             # operator for periodic boundary
    FETypes = [H1P1{size(xgrid[Coordinates],1)}],   # FETypes (default: P1)
    linsolver = ExtendableSparse.MKLPardiso,
    nsteps = ones(Int,length(FETypes)),             # number of embedding steps (parameters are scaled by nsteps equidistant steps within 0:1)
    damping = 0,                                    # damping in Newton iteration
    target_residual = 1e-12*ones(Int,length(FETypes)),
    solve_polarisation = false,                     # also solve for polarisation ?
    coupled = false,                                # solve coupled system in each nonlinear iteration ? (not needed as long as polarisation does not influence displacement)
    bonus_quadorder = 2,
    use_sparsity_detection = false,                 # sparsity detection for jacobians
    maxiterations = 20*ones(Int,length(FETypes)),
    fixed_penalty = 1e60
    )

    ## generate FES spaces
    Tv = Float64
    Ti = Int32
    FES = Array{FESpace{Tv,Ti},1}(undef, 0)
    push!(FES, FESpace{FETypes[1]}(xgrid))

    solve_polarisation = (length(FETypes) == 2) * solve_polarisation
    coupled = coupled * solve_polarisation

    if solve_polarisation
        push!(FES, FESpace{FETypes[2]}(xgrid))
    end

    ## generate FEVectors and FEMatrix
    npartitions = num_partitions(xgrid)
    SolutionOld = FEVector(FES)     # storage for previous solution
    Solution::FEVector{Float64} = FEVector(FES)        # storage for current and final solution
    if coupled
        S = FEMatrix([FES]; npartitions = npartitions)       # system matrix
        rhs = FEVector(FES)     # system rhs
    else
        S = FEMatrix([FES[1]]; npartitions = npartitions)       # system matrix
        rhs = FEVector(FES[1])     # system rhs
    end

    @info "Starting low level solver...
            ndofs = $(length(Solution.entries))
            polarisation = $solve_polarisation (coupled = $coupled)"

    ## prepare linear solver
	LP = LinearProblem(S.entries.cscmatrix, rhs.entries)
    linsolve = nothing

    ## prepare parallelization
    ncolors = num_pcolors(xgrid)
    pc = xgrid[PartitionCells]
    pcp = xgrid[PColorPartitions]
    itemassemblygroups = [pc[j]:(pc[j + 1] - 1) for j in 1:num_partitions(xgrid)]
    ngroups = length(itemassemblygroups)
    @info "ngroups = $ngroups, ncolors = $(num_pcolors(xgrid))"

    ## generate quadrature rule
    EG = xgrid[UniqueCellGeometries][1]
    order::Int = get_polynomialorder(FETypes[1], EG)
    qf = QuadratureRule{Float64, EG}(2*(order-1) + bonus_quadorder)
    weights::Vector{Float64} = qf.w
    nweights::Int = length(weights)
    cellvolumes::Array{Float64,1} = xgrid[CellVolumes]
    cellregions::Array{Int,1} = xgrid[CellRegions]

    ## prepare dofmaps and FE evaluators
    ndofs_u::Int = get_ndofs(ON_CELLS, FETypes[1], EG)      # displacement
    ∇u4p::Array{FEEvaluator,1} = [FEEvaluator(FES[1], Gradient, qf) for g=1:ngroups]
    celldofs_u::Matrix{Int32} = Matrix(FES[1][CellDofs])
    indices_∇u::Vector{Int} = 1:size(∇u4p[1].cvals, 1)

    if solve_polarisation # polarisation
        ndofs_P::Int = get_ndofs(ON_CELLS, FETypes[2], EG)  
        ∇P = [FEEvaluator(FES[2], Gradient, qf) for g=1:ngroups]
        vals_∇P = [∇P.cvals for g=1:ngroups]
        celldofs_P::Matrix{Int32} = Matrix(FES[2][CellDofs])
        indices_∇P::Vector{Int} = 1:size(vals_∇P, 1)
    end

    ## local matrix and vector structures
    AUU4p = [zeros(Float64, ndofs_u, ndofs_u) for g = 1:ngroups]
    bU4p = [zeros(Float64, ndofs_u) for g = 1:ngroups]

    SE = S.entries
    rhsE::Array{Float64,1} = rhs.entries
    SolutionE::Array{Float64,1} = Solution.entries
    offsets::Array{Int,1} = [0, FES[1].ndofs]
    if solve_polarisation
        push!(offsets, FES[1].ndofs + FES[2].ndofs)
    end

    ## prepare automatic differentation of displacement operator
    dim = size(xgrid[Coordinates],1)
    if solve_polarisation
        @assert dim == polarisation_operator.dim
    end
    argsizes_displacement::Vector{Int} = [dim^2,dim^2,Int(dim^2+dim*(dim+1)/2)]
    argsizes_polarisation::Vector{Int} = [dim,dim+dim^2,Int(dim+dim^2+dim*(dim+1)/2)]
    QP4p = [QPInfos(xgrid) for g = 1:ngroups]
    displacement_operator4qp4p = [(result,input) -> displacement_operator(result,input,QP) for QP in QP4p]
    polarisation_operator4qp4p = [(result,input) -> polarisation_operator(result,input,QP) for QP in QP4p]
    result_displacement4p = [Vector{Float64}(undef,argsizes_displacement[1]) for g = 1:ngroups]
    input_displacement4p = [Vector{Float64}(undef,argsizes_displacement[3]) for g = 1:ngroups]
    Dresult_displacement4p = [DiffResults.JacobianResult(result_displacement4p[g],input_displacement4p[g]) for g = 1:ngroups]
    result_polarisation4p = [Vector{Float64}(undef,argsizes_polarisation[1]) for g = 1:ngroups]
    input_polarisation4p = [Vector{Float64}(undef,argsizes_polarisation[3]) for g = 1:ngroups]
    Dresult_polarisation4p = [DiffResults.JacobianResult(result_polarisation4p[g],input_polarisation4p[g]) for g = 1:ngroups]

    if (use_sparsity_detection)
        sparsity_pattern = Symbolics.jacobian_sparsity(displacement_operator4qp,result_displacement,input_displacement)
        jac = Float64.(sparse(sparsity_pattern))
        colors = matrix_colors(jac) 
        cfg_displacement4p = [ForwardColorJacCache(displacement_operator4qp[g],input_displacement[g],nothing;
                dx = nothing,
                colorvec = colors,
                sparsity = nothing) for g = 1 :ngroups]
    else
        cfg_displacement4p = [ForwardDiff.JacobianConfig(displacement_operator4qp4p[1], result_displacement4p[g], input_displacement4p[g], ForwardDiff.Chunk{argsizes_displacement[3]}()) for g = 1:ngroups]
        cfg_polarisation = ForwardDiff.JacobianConfig(polarisation_operator4qp4p[1], result_polarisation4p[1], input_polarisation4p[1], ForwardDiff.Chunk{argsizes_polarisation[3]}())
    end
    
    ## ASSEMBLY LOOP
    bdofs = nothing
    bdofsP = nothing
    eval_∇u4p = [zeros(Float64, argsizes_displacement[3]) for g=1:ngroups]
    eval_∇P4p = [zeros(Float64, argsizes_polarisation[3]) for g=1:ngroups]
    ncells::Int = num_cells(xgrid)
    tempV4p = [zeros(Float64, argsizes_displacement[1]) for g=1:ngroups]
    tempP4p = [zeros(Float64, argsizes_polarisation[1]) for g=1:ngroups]


    function update_displacement_part(AUU, bU, QP, ∇u, Dresult_displacement, cfg_displacement, cells, part)
        eval_∇u::Vector{Float64} = eval_∇u4p[part]
        tempV::Vector{Float64} = tempV4p[part]
        vals_∇u::Array{Float64,3} = ∇u.cvals
        displacement_operator4qp = displacement_operator4qp4p[part]
        jac_displacement::Matrix{Float64} = DiffResults.jacobian(Dresult_displacement)
        jac_displacement_view = view(jac_displacement,indices_∇u,indices_∇u)
        value_displacement::Vector{Float64} = DiffResults.value(Dresult_displacement)

        for cell in cells

            ## update FE basis evaluators
            QP.item = cell
            QP.region = cellregions[cell]
            ∇u.citem[] = cell
            update_basis!(∇u) 
            dofs_u = view(celldofs_u, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u in current solution at current quadrature point
                fill!(eval_∇u, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇u[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_displacement, displacement_operator4qp, eval_∇u, cfg_displacement)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_displacement,displacement_operator4qp,value_displacement,eval_∇u,cfg_displacement)
                end

                # update matrix
                for j = 1 : ndofs_u
                    # multiply with jacobian
                    mul!(tempV, jac_displacement_view, view(vals_∇u,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_u
                        AUU[j,k] += dot(tempV, view(vals_∇u,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempV, jac_displacement_view, view(eval_∇u, indices_∇u))
                displacement_operator(value_displacement, eval_∇u, QP)
                tempV .-= value_displacement
                for j = 1 : ndofs_u
                    bU[j] += dot(tempV, view(vals_∇u,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_u
                for k = 1 : ndofs_u
                    if abs(AUU[j,k]) > 0
                        rawupdateindex!(SE, +, AUU[j,k] * cellvolumes[cell] , dofs_u[k], dofs_u[j], part)
                    end
                end
                rhsE[dofs_u[j]] += bU[j] * cellvolumes[cell]
            end

            fill!(AUU, 0)
            fill!(bU, 0)
        end
    end
    

    function update_displacement(; coupled = false)
        fill!(S[1,1],0)
        fill!(rhs[1],0)
        if coupled
            fill!(S[1,2],0)
        end
        
        for color in 1:ncolors
            Threads.@threads for part in pcp[color]:(pcp[color + 1] - 1)
                if (use_sparsity_detection)
                    sparsity_pattern = Symbolics.jacobian_sparsity(displacement_operator4qp,result_displacement,input_displacement)
                    jac = Float64.(sparse(sparsity_pattern))
                    cfg_displacement = ForwardColorJacCache(displacement_operator4qp[part],input_displacement[part],nothing;
                            dx = nothing,
                            colorvec = matrix_colors(jac) ,
                            sparsity = nothing)
                else
                    cfg_displacement = ForwardDiff.JacobianConfig(displacement_operator4qp4p[1], result_displacement4p[part], input_displacement4p[part], ForwardDiff.Chunk{argsizes_displacement[3]}())
                end
                update_displacement_part(AUU4p[part], bU4p[part], QP4p[part], ∇u4p[part], Dresult_displacement4p[part], cfg_displacement, itemassemblygroups[part], part)
            end
        end
        flush!(SE)

        ## apply penalties for boundary condition
        if boundary_operator !== nothing
            assemble!(boundary_operator, FES[1])
            bdofs = fixed_dofs(boundary_operator)
            for dof in bdofs
                rawupdateindex!(SE, +, fixed_penalty , dof, dof, 1)
                rhs.entries[dof] = 0
                Solution.entries[dof] = 0
            end
        end


        ## apply changes for periodic boundary
        if periodic_boundary_operator !== nothing
            ExtendableFEM.build_assembler!(periodic_boundary_operator, [Solution[1], Solution[1]])
            periodic_boundary_operator.assembler(SE, rhs.entries, true, true)
            bdofs = fixed_dofs(periodic_boundary_operator)
            for dof in bdofs
                SE[dof,dof] = fixed_penalty
                rhs.entries[dof] = 0
                Solution.entries[dof] = 0
            end
        end

        flush!(SE)
    end

    if solve_polarisation
        jac_polarisation_view = view(jac_polarisation,indices_∇P,dim^2+1:dim^2+dim)
        APP::Array{Float64,2} = zeros(Float64, ndofs_P, ndofs_P)
        bP::Array{Float64,1} = zeros(Float64, ndofs_P)
    end

    function update_polarisation(; coupled = false)

        if coupled
            offset = FES[1].ndofs
            fill!(S[2,1],0)
            fill!(S[2,2],0)
            fill!(rhs[2],0)
        else
            offset = 0
            fill!(S[1,1],0)
            fill!(rhs[1],0)
        end


        for cell = 1 : ncells

            ## update FE basis evaluators
            item_information[1] = cell
            item_information[2] = cell
            item_information[3] = cellregions[cell]
            pola_op = polarisation_operator4region(item_information)
            ∇u.citem[] = cell
            update_basis!(∇u) 
            ∇P.citem[] = cell
            update_basis!(∇P) 
            dofs_u = view(celldofs_u, : , cell)
            dofs_P = view(celldofs_P, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u and ∇P in current solution at current quadrature point
                fill!(eval_∇P, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇P[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end
                for j = 1 : ndofs_P, k in indices_∇P
                    eval_∇P[9+k] += SolutionE[dofs_P[j]+offsets[2]] * vals_∇P[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_polarisation, pola_op, eval_∇P, cfg_polarisation)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_polarisation,pola_op,result_polarisation,eval_∇P,cfg_polarisation)
                end

                # update matrix
                for j = 1 : ndofs_P
                    # multiply with jacobian
                    mul!(tempP, jac_polarisation_view, view(vals_∇P,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_P
                        APP[j,k] += dot(tempP, view(vals_∇P,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempP, jac_polarisation_view, view(eval_∇P, indices_∇P))
                pola_op(value_polarisation, eval_∇P)
                tempP .-= value_polarisation
                for j = 1 : ndofs_P
                    bP[j] += dot(tempP, view(vals_∇P,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_P
                for k = 1 : ndofs_P
                    if abs(APP[j,k]) > 0
                        rawupdateindex!(SE, +, APP[j,k] * cellvolumes[cell] , dofs_P[k], dofs_P[j])
                    end
                end
                rhs.entries[dofs_P[j]] += bP[j] * cellvolumes[cell]
            end

            fill!(APP, 0)
            fill!(bP, 0)
        end

        ## apply penalties for boundary condition
        fixed_dofsP = [offset + 1]
        for dof in fixed_dofsP
            SE[dof,dof] = fixed_penalty
            rhs.entries[dof] = 0 # fixed_penalty * Solution.entries[dof]
            Solution.entries[dof] = 0
        end

        flush!(SE)
    end

    ## prepare other stuff for loop
    emb_params_target = deepcopy(emb_params)
    factorization = nothing
    residual = zeros(Float64, size(SE,1))
    nlres::Float64 = 0
    linres::Float64 = 0
    change::Float64 = 0
    time_assembly::Float64 = 0
    time_solver::Float64 = 0
    iteration::Int = 0

    for step = 1 : nsteps[1]
        emb_params .= nsteps[1] == 1 ? emb_params_target : (step-1)/(nsteps[1]-1) .* emb_params_target

        @info "Entering embedding step = $step/$(nsteps[1]), emb_params = $emb_params"
            
        iteration = 0
        time_solver = 0
        while iteration <= maxiterations[1]
            # update system
            time_assembly = @elapsed update_displacement(; coupled = coupled)
            if iteration == 0 && step == 1
                linsolve = init(LP, linsolver)
            end
        
            # compute nonlinear residual
            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            if boundary_operator !== nothing
                view(residual, bdofs) .= 0
            end

            nlres = norm(residual)
            @info "         ---> nonlinear iteration $iteration, linres = $linres, nlres = $nlres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
            if nlres < target_residual[1]
                @info "target reached, stopping..."
                break
            end

            # update factorization and solve
            if damping > 0
                SolutionOld.entries .= Solution.entries
            end

            time_solver = @elapsed begin
				linsolve.A = S.entries.cscmatrix
				linsolve.b = rhs.entries
				x = LinearSolve.solve!(linsolve)
                Solution.entries .= x.u
            end
            iteration += 1

            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            if boundary_operator !== nothing
                view(residual, bdofs) .= 0
            end
            linres = norm(residual)

            if damping > 0
                view(Solution[1]) .= (1 - damping) * Solution.entries + damping * SolutionOld.entries
                change = norm(SolutionOld.entries - Solution.entries)
            end
        end
    end

    if !coupled && solve_polarisation
        @info "Entering decoupled polarisation solve..."
        S = FEMatrix(FES[2])       # system matrix
        rhs = FEVector(FES[2])     # system rhs
        SE = S.entries
        residual = zeros(Float64, size(SE,1))

        time_assembly = @elapsed update_polarisation(; coupled = coupled)
        factorization=linsolver(SE)
        factorize!(factorization, SE)

        time_solver = @elapsed begin
            ExtendableSparse.update!(factorization)
            ldiv!(view(Solution[2]), factorization, rhs.entries)
        end
        iteration += 1

        fill!(residual, 0)
        mul!(residual, SE, view(Solution[2]))
        residual .-= rhs.entries
        view(residual, fixed_dofsP) .= 0
        linres = norm(residual)

        @info "         ---> linres = $linres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
        
    end

    return Solution, nlres

end