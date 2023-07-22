module HydroTools

include("ET_models.jl")
include("heat_index.jl")
include("HW_index.jl")
include("detect_events.jl")


export cal_es, Tdew2RH, Tdew2VPD
export cal_U2, cal_lambda, cal_slope,
  ET0_eq, ET0_Penman48, ET0_FAO98, 
  heat_index


end # module HydroTools
