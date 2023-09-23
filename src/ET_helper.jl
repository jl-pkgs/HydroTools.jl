#' @param z.wind Height where measure the wind speed `[m]`. Default 10m.
#' @param Uz wind speed at height `z.wind`
#' @param U2 wind speed at 2m

"""

# Returns
- `cal_lambda` : "MJ / kg"
- `cal_slope`  : "kPa / K"
- `cao_rho_a`  : air density, "kg/m^3"
"""
function cal_U2(Uz::T, z=10.0) where {T<:Real}
  Uz * 4.87 / log(67.8 * z - 5.42)
end

function cal_Uz(U2::T, z)::T  where {T<:Real}
  log(67.8 * Z - 5.42) / 4.87 * U2
end


function cal_lambda(Tair::T) where {T<:Real}
  (2500.0 - Tair * 2.2) / 1000.0 * u"MJ / kg"
end

function cal_slope(Tair::T) where {T<:Real}
  4098.0 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) / (Tair + 237.3)^2 * u"kPa / K"
end


# rho_a: kg m-3
function cal_rho_a(Tair, q) 
  3.486 * Pa/cal_TvK(Tair, q) # FAO56, Eq. 3-5
end

## 几种计算虚温的方法

# FAO56, Eq. 3-7
cal_TvK(Tair) = 1.01 * (Tair + 273)

# 这个是最精确的版本
# FAO56, Eq. 3-6
cal_TvK(Tair, ea, Pa) = (Tair + K0) * (1 + (1 - epsilon) * ea / Pa)

# https://github.com/CUG-hydro/class2022_CUG_HydroMet/blob/master/ch02_大气的基本特征.md
# q ≈ ϵ*ea/Pa
# q = ϵ*ea/(Pa - (1 - ϵ)*ea)
function cal_TvK(Tair, q)
  # ea / Pa = q/(ϵ + (1 - ϵ) * q)
  (Tair + K0) * (1 + (1 - epsilon) * q / (ϵ + (1 - ϵ) * q))
end


W2mm(Ra, lambda) = Ra * 86400 / 1e6 / lambda

export cal_Uz, cal_U2, cal_lambda, cal_slope
