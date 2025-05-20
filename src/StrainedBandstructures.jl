module StrainedBandstructures

using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using SimplexGridFactory
using Triangulate
using TetGen
using LinearAlgebra
using LinearSolve
using Printf
using DiffResults
using ForwardDiff
using Symbolics
using SparseArrays
using SparseDiffTools
using DocStringExtensions



## problem: loading and saving grids leads to ElementGeometries -> DataType conversion (by DrWarson/JLD2?) which has to be reverted
## after loading (until this is fixed ina more elegant way)
function repair_grid!(xgrid::ExtendableGrid)
    xgrid.components[CellGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[CellGeometries][1], num_cells(xgrid))
    xgrid.components[FaceGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[FaceGeometries][1], length(xgrid.components[FaceGeometries]))
    xgrid.components[BFaceGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[BFaceGeometries][1], length(xgrid.components[BFaceGeometries]))

    xgrid.components[UniqueCellGeometries] = Vector{ElementGeometries}([xgrid.components[CellGeometries][1]])
    xgrid.components[UniqueFaceGeometries] = Vector{ElementGeometries}([xgrid.components[FaceGeometries][1]])
    xgrid.components[UniqueBFaceGeometries] = Vector{ElementGeometries}([xgrid.components[BFaceGeometries][1]])
end
export repair_grid!

include("grids.jl")
export bimetal_strip3D, bimetal_strip2D, bimetal_strip3D_middle_layer, bimetal_tensorgrid, bimetal_tensorgrid_uniform, bimetal_tensorgrid_uniform!
export condensator3D, condensator3D_tensorgrid, condensator2D
export nonpolarquantumwell3D, nonpolarquantumwell2D
export nanowire_grid, nanowire_tensorgrid, nanowire_tensorgrid_mirror

include("materialstructuretype.jl")
export MaterialStructureType
export ZincBlende2D, ZincBlende001, ZincBlende111_C14, ZincBlende111_C15, ZincBlende111_C14_C15, Wurtzite0001

include("materials.jl")
export MaterialType
export TestMaterial, GaAs, AlInAs, AlGaAs, InGaAs
export MaterialParameters
export get_materialtype

include("tensors.jl")
export AbstractTensor
export DenseMatrixTensor, SparseMatrixTensor
export IsotropicElasticityTensor
export apply_tensor!


include("heterostructures.jl")
export HeteroStructureData
export set_data

include("strains.jl")
export StrainType
export LinearStrain, LinearStrain2D, LinearStrain3D
export NonlinearStrain, NonlinearStrain2D, NonlinearStrain3D
export PreStrainType
export IsotropicPrestrain, AnisotropicDiagonalPrestrain
export eval_strain!, eval_elastic_strain!

include("pdeoperators.jl")
export PDEDisplacementOperator, PDEPolarisationOperator
export get_displacement_operator, get_displacement_operator
export get_polarisation_from_strain_operator
export get_polarisation_laplacian_operator
export get_polarisation_operator # wip
export get_energy_integrator

include("plane_cuts.jl")
export tet_x_plane!, plane_cut, get_cutgrids
export perform_simple_plane_cuts, perform_plane_cuts

include("export.jl")
export exportVTK

include("statistics.jl")
export compute_statistics

include("solvers.jl")
export solve_by_embedding
export solve_by_embedding!
export solve_by_damping
export solve_lowlevel, solve_lowlevel_parallel

include("solvers_from_energy.jl")
export EnergyOperator, update_M!
export ∂FW!

greet() = print("Hello World!")

end # module
