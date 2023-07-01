include("constant.jl")
include("humidity.jl")

#' @param z.wind Height where measure the wind speed `[m]`. Default 10m.
#' @param Uz wind speed at height `z.wind`
#' @param U2 wind speed at 2m
function cal_U2(Uz::T, z_wind=10.0)::T where {T<:Real}
  z_wind == 2.0 ? Uz : Uz * 4.87 / log(67.8 * z_wind - 5.42)
end

function cal_lambda(Tair::T)::T where {T<:Real} 
  (2500.0 - Tair * 2.2) / 1000.0
end

function cal_slope(Tair::T)::T where {T<:Real}
  4098.0 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) / (Tair + 237.3)^2
end

"""

# Arguments
- `Rn`   : net radiation (W m-2)
- `Tair` : 2m air temperature (degC)
- `VPD`  : vapor pressure deficit (kPa)
- `Pa`   : surface air pressure (kPa)
"""
function ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...) where {T<:Real}

  lambda::T = cal_lambda(Tair) # MJ kg-1
  slope::T = cal_slope(Tair) # kPa degC-1
  gamma::T = Cp * Pa / (epsilon * lambda) # kPa degC-1

  coef_W2mm::T = 0.086400 / lambda
  Eeq::T = slope / (slope + gamma) * Rn * coef_W2mm

  lambda, slope, gamma, Eeq
end


# ET0_FAO98(200.0, 20.0, 2.0, 2.0)
function ET0_Penman48(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2) where {T<:Real}
  lambda, slope, gamma, Eeq = ET0_eq(Rn, Tair, Pa)

  U2::T = cal_U2(wind, z_wind)
  # rho_a * Cp * dT / rH (MJ m-2 s-1)
  # rho_a ≈ 1.225 kg/m3
  # rho_a * Cp / rH = f(U2)
  # `f(U2) = 2.6 * (1 + 0.54U2)` is equivalent to Shuttleworth1993
  Evp::T = gamma / (slope + gamma) * 6.43 * (1 + 0.536 * U2) * VPD / lambda
  PET::T = Eeq + Evp
  PET
end


# ET0_Penman48(200., 20., 2., 2.)
function ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false) where {T<:Real}
  U2 = cal_U2(wind, z_wind)
  
  if tall_crop
    p1 = T(1600.0)
    p2 = T(0.38)
  else
    p1 = T(900.0)
    p2 = T(0.34)
  end

  lambda, slope, gamma, Eeq = ET0_eq(Rn, Tair, Pa)
  coef_W2mm::T = 0.086400 / lambda

  Eeq::T = slope / (slope + (gamma * (1.0 + p2 * U2))) * Rn * coef_W2mm
  Evp::T = gamma * p1 / (Tair + 273.15) * U2 * VPD / (slope + (gamma * (1.0 + p2 * U2)))

  PET::T = Eeq + Evp
  PET
end
