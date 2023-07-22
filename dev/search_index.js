var documenterSearchIndex = {"docs":
[{"location":"ExtremeClimate/#Extreme-Climate-indexes","page":"Extreme Climate indexes","title":"Extreme Climate indexes","text":"","category":"section"},{"location":"ExtremeClimate/#Events-detection","page":"Extreme Climate indexes","title":"Events detection","text":"","category":"section"},{"location":"ExtremeClimate/","page":"Extreme Climate indexes","title":"Extreme Climate indexes","text":"extract_eventInfo\ndetect_events","category":"page"},{"location":"ExtremeClimate/#HydroTools.extract_eventInfo","page":"Extreme Climate indexes","title":"HydroTools.extract_eventInfo","text":"extract_eventInfo(y::AbstractVector{<:Real}, i_beg::Integer, i_end::Integer;\n    index=nothing, goal=1, len_min=1, len2peak=0, ignored...)\n\nExtracts information about an event in a signal.\n\nExamples\n\njulia> y = [0, 1, 2, 3, 2, 1, 0];\njulia> extract_eventInfo(y, 2, 6)\n(i_beg = 2, i_peak = 4, i_end = 5, len = 4, len_left = 2, len_right = 1, peak = 3)\n\n\n\n\n\n","category":"function"},{"location":"ExtremeClimate/#HydroTools.detect_events","page":"Extreme Climate indexes","title":"HydroTools.detect_events","text":"detect_events(lgl::AbstractVector{Bool}; len_min=1, index=nothing, only_max=false, ignored...)\ndetect_events(y::AbstractVector{<:Real}, lgl::BitVector;\n    len_min=1,\n    len2peak=0,\n    goal=1,\n    index=nothing,\n    only_max=false,\n    ignored...\n)\n\nDetects events in a signal based on a logical vector.\n\nArguments\n\ny::AbstractVector{<:Real}: The signal to detect events in.\nlgl::BitVector: A logical vector indicating where events occur in the signal.\nlen_min=1: (optional) The minimum length of an event.\nlen2peak=0: (optional) The minimum length from the peak to the beginning or end of an event.\ngoal=1: (optional), -1 or 1. If 1, find the maximum value in the event; If -1, minimum value used.\nindex=nothing: (optional) The index of the signal. If provided, the returned indices will be in the index space.\nonly_max=false: (optional) If true, only events with the maximum duration will be returned.\nignored...: (optional) Ignored arguments.\n\nReturns\n\nAn array of named tuples, where each tuple represents an event and has the following fields:\n\ni_beg::Integer: The index of the beginning of the event.\ni_peak::Integer: The index of the peak of the event.\ni_end::Integer: The index of the end of the event.\nlen::Integer: The length of the event.\nlen_left::Integer: The length from the beginning of the event to the peak.\nlen_right::Integer: The length from the peak to the end of the event.\npeak::Real: The value of the peak of the event.\n\nIf only_max is true, only the event with the maximum duration will be returned.\n\nExamples\n\njulia> y = [0, 1, 2, 3, 0, 2, 1, 0];\njulia> lgl = y .> 0;\njulia> detect_events(lgl)\n2-element Vector{NamedTuple{(:i_beg, :i_end, :len), Tuple{Int64, Int64, Int64}}}:\n (i_beg = 2, i_end = 4, len = 3)\n (i_beg = 6, i_end = 7, len = 2)\njulia> detect_events(y, lgl)\n2-element Vector{NamedTuple{(:i_beg, :i_peak, :i_end, :len, :len_left, :len_right, :peak), NTuple{7, Int64}}}:\n (i_beg = 2, i_peak = 4, i_end = 4, len = 3, len_left = 2, len_right = 0, peak = 3)\n (i_beg = 6, i_peak = 6, i_end = 7, len = 2, len_left = 0, len_right = 1, peak = 2)\n\n\n\n\n\n","category":"function"},{"location":"ExtremeClimate/#Heatwave","page":"Extreme Climate indexes","title":"Heatwave","text":"","category":"section"},{"location":"ExtremeClimate/","page":"Extreme Climate indexes","title":"Extreme Climate indexes","text":"(Image: The definition of Heatwave events)","category":"page"},{"location":"ExtremeClimate/","page":"Extreme Climate indexes","title":"Extreme Climate indexes","text":"HW_index","category":"page"},{"location":"ExtremeClimate/#HydroTools.HW_index","page":"Extreme Climate indexes","title":"HydroTools.HW_index","text":"HW_index(anorm::AbstractVector; p_left = 0.99)\n\nCompute the HW index for a given anomaly vector anorm.\n\nArguments\n\nanorm::AbstractVector: A vector of anomaly scores.\np_left::Float64=0.99: The probability of false positives.\n\nReturns\n\nA named tuple with the following fields:\n\nduration::Int: The duration of the anomaly.\nfrequency::Int: The number of anomaly events.\nintensity::Float64: The maximum anomaly score.\nvolume::Float64: The sum of anomaly scores.\nPR::Float64: The probability of detection.\nFAR::Float64: The false alarm rate.\n\nExample\n\njulia> anorm = [0.1, 0.2, 0.3, 0.2, 0.1, 0, -0.1, 0.1, 0.2, 0.3]\njulia> HW_index(anorm)\n(duration = 9, frequency = 2, intensity = 0.3, volume = 1.5, PR = 89.99999999999993, FAR = 0.9888888888888889)\n\n\n\n\n\n","category":"function"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"] add https://github.com/CUG-hydro/HydroTools.jl","category":"page"},{"location":"#Contents","page":"Introduction","title":"Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Meteorological basics  \nHumidity\nRadiation\nPotential Evapotranspiration models\nExtreme Climate indexes\nHeatwave index","category":"page"},{"location":"#References","page":"Introduction","title":"References","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"https://github.com/CUG-hydro/hydroTools","category":"page"},{"location":"PET/#Potential-Evapotranspiration-models","page":"Potential Evapotranspiration models","title":"Potential Evapotranspiration models","text":"","category":"section"},{"location":"PET/","page":"Potential Evapotranspiration models","title":"Potential Evapotranspiration models","text":"ET0_eq\nET0_Penman48\nET0_FAO98","category":"page"},{"location":"PET/#HydroTools.ET0_eq","page":"Potential Evapotranspiration models","title":"HydroTools.ET0_eq","text":"ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...)\n\nEquilibrium evaporation, ET0_eq = slope / (slope + gamma) * Rn\n\nArguments\n\nRn   : net radiation (W m-2)\nTair : 2m air temperature (degC)\nVPD  : vapor pressure deficit (kPa)\nPa   : surface air pressure (kPa)\n\nReturns\n\nlambda : latent heat of vaporization (MJ kg-1)\nslope  : slope of the saturation vapor pressure curve (kPa degC-1)\ngamma  : psychrometric constant (kPa degC-1)\nEeq    : equilibrium evaporation rate (mm day-1)\n\nOptional keyword arguments\n\nargs... : additional arguments (not used in this function)\n\nExamples\n\njulia> ET0_eq(200.0, 20.0, 2.0)\n(2.456, 0.14474018811241365, 0.0013262323104019538, 6.971947723218883)\n\n\n\n\n\n","category":"function"},{"location":"PET/#HydroTools.ET0_Penman48","page":"Potential Evapotranspiration models","title":"HydroTools.ET0_Penman48","text":"ET0_Penman48(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2)\n\nExamples\n\nET0_Penman48(200., 20., 2., 2.)\n\n\n\n\n\n","category":"function"},{"location":"PET/#HydroTools.ET0_FAO98","page":"Potential Evapotranspiration models","title":"HydroTools.ET0_FAO98","text":"ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false)\n\nExamples\n\nET0_FAO98(200.0, 20.0, 2.0, 2.0)\n\n\n\n\n\n","category":"function"}]
}
