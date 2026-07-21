using NetCDFTools

uv2wind(u::T, v::T) where {T<:Real} = sqrt(u^2 + v^2)
convertU2(u10::T, v10::T) where {T<:Real} = T(cal_U2(uv2wind(u10, v10)))

"""
    _FAO98(Rn_G, Ta, VPD, U2, Pa=atm)

Hourly FAO-56 reference evapotranspiration (mm h⁻¹). `Rn_G` is hourly net
radiation minus soil heat flux (MJ m⁻² h⁻¹).
"""
function _FAO98(Rn_G::T, Ta::T, VPD::T, U2::T, Pa::T=T(atm)) where {T<:AbstractFloat}
  γ::T = cal_gamma(Ta, Pa)
  Δ::T = cal_slope(Ta)
  λ::T = cal_lambda(Ta)
  den::T = Δ + γ * (one(T) + T(0.34) * U2)

  PET_e::T = Δ / den * Rn_G / λ
  PET_a::T = γ * T(37) / (Ta + T(273.15)) * U2 * max(VPD, zero(T)) / den
  max(PET_e + PET_a, zero(T))
end

"""
    _PT1972(Rn_G, Ta, Pa=atm)

Hourly Priestley–Taylor evapotranspiration (mm h⁻¹). `Rn_G` is in
MJ m⁻² h⁻¹.
"""
function _PT1972(Rn_G::T, Ta::T, Pa::T=T(atm)) where {T<:AbstractFloat}
  γ::T = cal_gamma(Ta, Pa)
  Δ::T = cal_slope(Ta)
  λ::T = cal_lambda(Ta)
  PET::T = T(1.26) * Δ / (Δ + γ) * Rn_G / λ
  max(PET, zero(T))
end

"""
    _Romanenko(Ta, Tdew)

Approximate hourly Romanenko evapotranspiration (mm h⁻¹), obtained by dividing
the empirical daily rate by 24. Prefer daily-mean forcing for daily estimates.
"""
function _Romanenko(Ta::T, Tdew::T) where {T<:AbstractFloat}
  ea::T = cal_es(Tdew)
  es::T = cal_es(Ta)
  PET_day::T = T(4.5) * (one(T) + Ta / T(25))^2 * (one(T) - ea / es)
  max(PET_day / T(24), zero(T))
end

function ET0_FAO98(f::String)
  Tdew = nc_read(f, "dewpoint_temperature_2m") .- 273.15f0
  Ta = nc_read(f, "temperature_2m") .- 273.15f0
  Pa = nc_read(f, "surface_pressure") ./ 1000
  H = nc_read(f, "surface_sensible_heat_flux") ./ Float32(1e6)
  LE = nc_read(f, "surface_latent_heat_flux") ./ Float32(1e6)
  Rn_G = -(H + LE)
  clamp!(Rn_G, 0, Inf)

  u10 = nc_read(f, "u_component_of_wind_10m")
  v10 = nc_read(f, "v_component_of_wind_10m")
  U2 = convertU2.(u10, v10)
  VPD = cal_es.(Ta) - cal_es.(Tdew)

  _FAO98.(Rn_G, Ta, VPD, U2, Pa)
end

function ET0_PriestleyTaylor1972(f::String)
  Ta = nc_read(f, "temperature_2m") .- 273.15f0
  Pa = nc_read(f, "surface_pressure") ./ 1000
  H = nc_read(f, "surface_sensible_heat_flux") ./ Float32(1e6)
  LE = nc_read(f, "surface_latent_heat_flux") ./ Float32(1e6)
  Rn_G = -(H + LE)

  _PT1972.(Rn_G, Ta, Pa)
end

function ET0_Romanenko(f::String)
  Ta = nc_read(f, "temperature_2m") .- 273.15f0
  Tdew = nc_read(f, "dewpoint_temperature_2m") .- 273.15f0
  _Romanenko.(Ta, Tdew)
end

export ET0_FAO98, ET0_PriestleyTaylor1972, ET0_Romanenko
