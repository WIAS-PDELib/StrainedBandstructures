
"""
    $(TYPEDEF)

structure that collects material data sets, elasticity tensors and piezoelectric tensors for all materials
"""
struct HeteroStructureData{T, MST <: MaterialStructureType}
    data::Array{MaterialParameters,1}
    TensorC::Array{AbstractTensor{T},1}         # elastitity tensor C
    TensorE::Array{AbstractTensor{T},1}   # piezo-electricity tensor E
end

"""
````
function HeteroStructureData(materials::Array{MaterialParameters,1}, MST::Type{<:MaterialStructureType})
````

construct a HeteroStructureData from the given array of MaterialParameters and
generates the elasticity and piezo-eletricity tensors for the given MaterialStructureTypes
"""
function HeteroStructureData(materials::Array{DataType,1}, MST::Type{<:MaterialStructureType}, TensorMatrixType = DenseMatrixTensor)

    data = Array{MaterialParameters,1}(undef, length(materials))
    C = Array{AbstractTensor{Float64},1}(undef,length(materials))
    E = Array{AbstractTensor{Float64},1}(undef,length(materials))

    @info "using TensorMatrixType = $TensorMatrixType"
    for m = 1 : length(materials)

        ## get material data set
        data[m] = MaterialParameters(materials[m], MST)

        ## setup tensors
        if MST <: ZincBlende001
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]

            # stiffness tensor in N/(nm)^2 # todo: move scaling to solver
            # Equation (4)
            # in "Symmetry-adapted calculations of strain and polarization fields in (111)-oriented zinc-blende quantum dots" by Schulz et. al
            C[m] = TensorMatrixType(
              1e-9*[ C11 C12 C12   0   0   0
                     C12 C11 C12   0   0   0
                     C12 C12 C11   0   0   0
                       0   0   0 C44   0   0
                       0   0   0   0 C44   0
                       0   0   0   0   0 C44 ])

            # piezoelectric tensor in C/(nm^2)
            # Equation (22)
            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E[m] = TensorMatrixType(
                  1e-9*[0 0 0 E14zb     0     0
                   0 0 0     0 E14zb     0
                   0 0 0     0     0 E14zb])
        elseif MST <: ZincBlende2D
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]

            # stiffness tensor in N/(nm)^2
            C[m] = TensorMatrixType(
                   1e-9*[ C11 C12   0
                          C12 C11   0
                            0   0 C44 ])

            # piezoelectric tensor in C/(nm^2)
            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E[m] = TensorMatrixType(
                1e-9*[E14zb     0 0 0 0 0
                     0 E14zb 0 0 0 0])
        elseif MST <: ZincBlende111_C14 || MST <: ZincBlende111_C15 || MST <: ZincBlende111_C14_C15
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            sr2 = sqrt(2)
            C11p = (1/2)*(C11 + C12) + C44
            C12p = (1/6)*(C11 + 5*C12) - (1/3)*C44
            C44p = (1/3)*(C11 - C12 + C44)
            C33p = (3/2)*C11p - (1/2)*C12p - C44p
            C13p = -(1/2)*C11p + (3/2)*C12p + C44p
            C14p = sr2/6*(-C11 + C12 + 2*C44)
            C15p = (1/sr2)*C11p - (1/sr2)*C12p - sr2*C44p
            C66p = (1/2)*(C11p - C12p)

            # stiffness tensor in N/(nm)^2
            # Equation (14)
            if MST <: ZincBlende111_C14
                C[m] = TensorMatrixType(
                    1e-9*[  C11p    C12p    C13p    C14p    0       0
                            C12p    C11p    C13p    -C14p   0       0
                            C13p    C13p    C33p    0       0       0
                            C14p    -C14p   0       C44p    0       0
                            0       0       0       0       C44p    C14p
                            0       0       0       0       C14p    C66p])
            elseif MST <: ZincBlende111_C15
                C[m] = TensorMatrixType(
                    1e-9*[  C11p    C12p    C13p    0       C15p    0
                            C12p    C11p    C13p    0       -C15p   0
                            C13p    C13p    C33p    0       0       0
                            0       0       0       C44p    0       -C15p
                            C15p    -C15p   0       0       C44p     0
                            0       0       0       -C15p   0       C66p])
            elseif MST <: ZincBlende111_C14_C15
                C[m] = TensorMatrixType(
                    1e-9*[  C11p    C12p    C13p    C14p    C15p    0
                            C12p    C11p    C13p    -C14p   -C15p   0
                            C13p    C13p    C33p    0       0       0
                            C14p    -C14p   0       C44p    0       -C15p
                            C15p    -C15p   0       0       C44p     C14p
                            0       0       0       -C15p   C14p    C66p])
            end

            # piezoelectric tensor in C/(nm^2)
            # Equation (27)
            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E11 = - sqrt(2/3) * E14zb
            E12 = sqrt(2/3) * E14zb
            E15 = - sqrt(1/3) * E14zb
            E31 = E15
            E33 = 2/sqrt(3) * E14zb

            E[m] = TensorMatrixType(
                   1e-9*[E11 E12 0   0   E15 0
                    0   0   0   E15 0   E12
                    E31 E31 E33 0   0   0])
        elseif MST <: Wurtzite0001
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C11wz = (1/6) * (3*C11 + 3 * C12 + 6 * C44)
            C12wz = (1/6) * (1*C11 + 5 * C12 - 2 * C44)
            C13wz = (1/6) * (2*C11 + 4 * C12 - 4 * C44)
            C33wz = (1/6) * (2*C11 + 4 * C12 + 8 * C44)
            C44wz = (1/6) * (2*C11 - 2 * C12 + 2 * C44)
            C66wz = (1/6) * (1*C11 - 1 * C12 + 4 * C44)

            # stiffness tensor in N/(nm)^2
            # Equation (7)
            C[m] = TensorMatrixType(
                   1e-9*[ C11wz C12wz C13wz 0     0     0
                          C12wz C11wz C13wz 0     0     0
                          C13wz C13wz C33wz 0     0     0
                          0     0     0     C44wz 0     0
                          0     0     0     0     C44wz 0
                          0     0     0     0     0     C66wz ])

            # piezoelectric tensor in C/(nm^2)
            # Equation (25)
            E31wz = data[m].PiezoElectricConstants["E31wz"]
            E33wz = data[m].PiezoElectricConstants["E33wz"]
            E15wz = data[m].PiezoElectricConstants["E15wz"]
            E[m] = TensorMatrixType(
                   1e-9*[     0     0     0     0  E15wz  0
                         0     0     0 E15wz      0  0
                     E31wz E31wz E33wz     0      0  0 ])
                elseif MST <: Wurtzite0001
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C11wz = (1/6) * (3*C11 + 3 * C12 + 6 * C44)
            C12wz = (1/6) * (1*C11 + 5 * C12 - 2 * C44)
            C13wz = (1/6) * (2*C11 + 4 * C12 - 4 * C44)
            C33wz = (1/6) * (2*C11 + 4 * C12 + 8 * C44)
            C44wz = (1/6) * (2*C11 - 2 * C12 + 2 * C44)
            C66wz = (1/6) * (1*C11 - 1 * C12 + 4 * C44)

            # stiffness tensor in N/(nm)^2
            # Equation (7)
            C[m] = TensorMatrixType(
                   1e-9*[ C11wz C12wz C13wz 0     0     0
                          C12wz C11wz C13wz 0     0     0
                          C13wz C13wz C33wz 0     0     0
                          0     0     0     C44wz 0     0
                          0     0     0     0     C44wz 0
                          0     0     0     0     0     C66wz ])

            # piezoelectric tensor in C/(nm^2)
            # Equation (25)
            E31wz = data[m].PiezoElectricConstants["E31wz"]
            E33wz = data[m].PiezoElectricConstants["E33wz"]
            E15wz = data[m].PiezoElectricConstants["E15wz"]
            E[m] = TensorMatrixType(
                   1e-9*[     0     0     0     0  E15wz  0
                         0     0     0 E15wz      0  0
                     E31wz E31wz E33wz     0      0  0 ])
        elseif MST <: Wurtzite0001
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C11wz = (1/6) * (3*C11 + 3 * C12 + 6 * C44)
            C12wz = (1/6) * (1*C11 + 5 * C12 - 2 * C44)
            C13wz = (1/6) * (2*C11 + 4 * C12 - 4 * C44)
            C33wz = (1/6) * (2*C11 + 4 * C12 + 8 * C44)
            C44wz = (1/6) * (2*C11 - 2 * C12 + 2 * C44)
            C66wz = (1/6) * (1*C11 - 1 * C12 + 4 * C44)

            # stiffness tensor in N/(nm)^2
            # Equation (7)
            C[m] = TensorMatrixType(
                   1e-9*[ C11wz C12wz C13wz 0     0     0
                          C12wz C11wz C13wz 0     0     0
                          C13wz C13wz C33wz 0     0     0
                          0     0     0     C44wz 0     0
                          0     0     0     0     C44wz 0
                          0     0     0     0     0     C66wz ])

            # piezoelectric tensor in C/(nm^2)
            # Equation (25)
            E31wz = data[m].PiezoElectricConstants["E31wz"]
            E33wz = data[m].PiezoElectricConstants["E33wz"]
            E15wz = data[m].PiezoElectricConstants["E15wz"]
            E[m] = TensorMatrixType(
                   1e-9*[     0     0     0     0  E15wz  0 
                         0     0     0 E15wz      0  0 
                     E31wz E31wz E33wz     0      0  0 ])
            elseif MST <: Wurtzite
                C11wz = data[m].ElasticConstants["C11"]
                C12wz = data[m].ElasticConstants["C12"]
                C13wz = data[m].ElasticConstants["C13"]
                C33wz = data[m].ElasticConstants["C33"]
                C44wz = data[m].ElasticConstants["C44"]
                C66wz = (1/2) * (C11wz - C12wz)

                # stiffness tensor in N/(nm)^2
                # Equation (7)
                C[m] = TensorMatrixType(
                        1e-9*[ C11wz C12wz C13wz 0     0     0
                                C12wz C11wz C13wz 0     0     0
                                C13wz C13wz C33wz 0     0     0
                                0     0     0     C44wz 0     0
                                0     0     0     0     C44wz 0
                                0     0     0     0     0     C66wz ])

                # piezoelectric tensor in C/(nm^2)
                # Equation (25)
                E31wz = data[m].PiezoElectricConstants["E31wz"]
                E33wz = data[m].PiezoElectricConstants["E33wz"]
                E15wz = data[m].PiezoElectricConstants["E15wz"]
                E[m] = TensorMatrixType(
                        1e-9*[     0     0     0     0  E15wz  0
                                0     0     0 E15wz      0  0
                            E31wz E31wz E33wz     0      0  0 ])
        end
    end

    return HeteroStructureData{Float64,MST}(data,C,E)
end
