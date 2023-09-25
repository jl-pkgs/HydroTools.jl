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


function cal_lambda(Tair::T) where {T<:Real}
  #  * u"MJ / kg"
  (2500.0 - Tair * 2.2) / 1000.0
end

function cal_slope(Tair::T) where {T<:Real}
  #  * u"kPa / K"
  4098.0 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) / (Tair + 237.3)^2
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


"""
    aerodynamic_conductance(U2, hc)

# Arguments
- `U2`: wind speed at 2m
- `hc`: canopy height

# Return
- `Ga`: aerodynamic conductance in m/s
"""
function aerodynamic_conductance(U2::T, hc::T; Zob = 15.0) where {T<:Real}
  kmar = 0.40        # von Karman's constant 0.40
  d = 0.64 * hc
  zom = 0.13 * hc
  zoh = 0.10 * zom
  uz = cal_Uz(U2, Zob) # convert from u2 to uz
  @fastmath Ga = uz * kmar^2 / (log((Zob - d) / zom) * log((Zob - d) / zoh)) # m s-1
  Ga
end


export cal_Uz, cal_U2, cal_lambda, cal_slope, mol2m
