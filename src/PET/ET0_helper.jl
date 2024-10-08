"""

# Returns
- `cal_lambda` : "MJ / kg"
- `cal_slope`  : "kPa / K"
- `cao_rho_a`  : air density, "kg/m^3"
"""
function cal_U2(Uz::T, z=10.0) where {T<:Real}
  z == 2 && (return Uz)
  Uz * 4.87 / log(67.8 * z - 5.42)
end

function cal_Uz(U2::T, z)::T where {T<:Real}
  z == 2 && (return U2)
  log(67.8 * z - 5.42) / 4.87 * U2
end

cal_Pa(z::T) where {T<:Real} = 
  atm * (((293 - 0.0065 * z) / 293)^5.26) # kPa, Allen 1998 Eq. 7

function cal_lambda(Tair::T) where {T<:Real}
  #  * u"MJ / kg"
  2.501 - 0.00237 * Tair # bolton 1980
  # (2500.0 - Tair * 2.2) / 1000.0
end

function cal_slope(Tair::T) where {T<:Real}
  #  * u"kPa / K"
  4098.0 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) / (Tair + 237.3)^2
end

function cal_gamma(Tair::T, Pa=atm) where {T<:Real}
  #  * u"kPa / K"
  lambda = cal_lambda(Tair)
  Cp * Pa / (ϵ * lambda)
end

function cal_bowen(Tair::T, Pa=atm) where {T<:Real}
  Δ = cal_slope(Tair)
  γ = cal_gamma(Tair, Pa)
  γ / Δ
end

# rho_a: kg m-3
function cal_rho_a(Tair, q, Pa)
  3.486 * Pa / cal_TvK(Tair, q) # FAO56, Eq. 3-5
end

cal_rho_a(Tair, Pa) = 3.486 * Pa / cal_TvK(Tair)


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


export cal_Uz, cal_U2, cal_Pa, cal_lambda, cal_gamma, cal_slope, cal_bowen, mol2m
