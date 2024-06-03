# using ProtoStructs
include("Met.jl")
include("Radiation.jl")
include("Canopy.jl")
include("Flux.jl")

include("soil_thermal_properties.jl")
include("Soil.jl")

# 三个高度：ref, surf, ground：参考高度，0表面，
export Met, Radiation, Canopy, Soil, Flux
export cal_coszen
