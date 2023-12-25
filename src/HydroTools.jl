module HydroTools

using Dates
using Dates: year
using DataStructures: OrderedDict

include("GOF.jl")
include("Optim/Optim.jl")

include("constant.jl")

include("cal_sun_angle.jl")
include("cal_radiation.jl")
include("cal_humidity.jl")
include("thermal.jl")

include("ET0_helper.jl")
include("ET0_models.jl")


include("ET/ET.jl")

include("unit_convert.jl")
include("heat_index.jl")
include("HW_index.jl")
include("detect_events.jl")
include("Climate/ClimateIndex.jl")

include("Radiation/Radiation.jl")

export cal_es, Tdew2RH, Tdew2VPD
export cal_U2, cal_lambda, cal_slope,
  ET0_eq, ET0_Penman48, ET0_FAO98, 
  heat_index


end # module HydroTools
