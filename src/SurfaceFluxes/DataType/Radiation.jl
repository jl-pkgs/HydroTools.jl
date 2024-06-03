"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Radiation{T<:AbstractFloat}
  "Surface emissivity"
  ϵ::T = 0.98
  "fraction of snow"
  fsnow::T = 0.0
  "albedo of surface and snow for visible and near-infrared wavebands"
  α_surf_vis::T = 0.10
  α_surf_nir::T = 0.20
  α_snow_vis::T = 0.95
  α_snow_nir::T = 0.70
  α_vis::T = α_surf_vis * (1 - fsnow) + α_snow_vis * fsnow
  α_nir::T = α_surf_nir * (1 - fsnow) + α_snow_nir * fsnow

  "sin of solar zenith angle"
  coszen::T = T(NaN)

  "Longwave radiation from the atmosphere, [W m-2]"
  Rln_in::T = T(NaN)
  "Clear sky total solar radiation at the top of the atmosphere, [W m-2]"
  Rs_toa::T = T(NaN)
  "Clear sky direct beam, [W m-2]"
  Rs_dir::T = T(NaN)
  "Clear sky diffuse, [W m-2]"
  Rs_dif::T = T(NaN)
  "Clear sky total solar radiation, [W m-2]"
  Rs::T = Rs_dir + Rs_dif
  "Inward Radiation, [W m-2]"
  Qa::T = (1 - α_vis) * 0.5 * Rs + (1 - α_nir) * 0.5 * Rs + ϵ * Rln_in
  "Net radiation, [W m-2]"
  Rn::T = T(NaN)
end

"""
# References
1. Gates, D.M. (1980), Clear sky atmospheric attenuation, Biophysical Ecology,
   110-115
"""
function Radiation(doy, hour, lat; fsnow=0.0, kw...)
  coszen = cal_coszen(doy, hour, lat)
  Rs_toa = cal_Rs_toa(doy, coszen) # Rs_toa

  τ_atm = 0.5
  oam = 1 / max(coszen, 0.04)
  Rs_dir = Rs_toa * τ_atm^oam # Clear sky direct beam
  Rs_dif = Rs_toa * (0.271 - 0.294 * τ_atm^oam) # Clear sky diffuse
  Rs = Rs_dif + Rs_dir # Clear sky total

  Radiation{Float64}(; Rs_toa, Rs, Rs_dir, Rs_dif, coszen,
    fsnow, kw...)
end

function cal_coszen(doy, hour, lat)
  decl = 23.45 * sin((284 + doy) / 365 * 2 * pi) * pi / 180
  hour_angle = 15 * (hour - 12) * pi / 180
  max(cos(lat) * cos(decl) * cos(hour_angle) + sin(lat) * sin(decl), 0) # coszen
end

function cal_Rs_toa(doy, coszen)
  Rg = 1364 # Solar radiation at top of the atmosphere, W/m2
  rv = 1 / sqrt(1 + 0.033 * cos(doy / 365 * 2 * pi))
  Rg / rv^2 * coszen # Rs_toa
end

function cal_Rs_toa(doy, hour, lat)
  coszen = cal_coszen(doy, hour, lat)
  cal_Rs_toa(doy, coszen) # Rs_toa
end
