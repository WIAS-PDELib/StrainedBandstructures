
"""
    $(TYPEDEF)
"""
abstract type MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type ZincBlende2D <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type ZincBlende001 <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type ZincBlende111_C14 <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type ZincBlende111_C15 <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type ZincBlende111_C14_C15 <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type Wurtzite0001 <: MaterialStructureType end
"""
    $(TYPEDEF)
"""
abstract type Wurtzite <: MaterialStructureType end
Base.String(::Type{ZincBlende001}) = "ZincBlende001"
Base.String(::Type{ZincBlende111_C14}) = "ZincBlende111_C14"
Base.String(::Type{ZincBlende111_C15}) = "ZincBlende111_C14"
Base.String(::Type{ZincBlende111_C14_C15}) = "ZincBlende111_C14_C15"
Base.String(::Type{Wurtzite0001}) = "Wurtzite0001"
Base.String(::Type{Wurtzite}) = "Wurtzite"
