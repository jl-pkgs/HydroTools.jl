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
  coszen::T = NaN

  "Longwave radiation from the atmosphere, [W m-2]"
  Rln_in::T = NaN
  "Clear sky total solar radiation at the top of the atmosphere, [W m-2]"
  Rs_toa::T = NaN
  "Clear sky direct beam, [W m-2]"
  Rs_dir::T = NaN
  "Clear sky diffuse, [W m-2]"
  Rs_dif::T = NaN
  "Clear sky total solar radiation, [W m-2]"
  Rs::T = Rs_dir + Rs_dif
  "Inward Radiation, [W m-2]"
  Qa::T = (1 - α_vis) * 0.5 * Rs + (1 - α_nir) * 0.5 * Rs + ϵ * Rln_in
  "Net radiation, [W m-2]"
  Rn::T = NaN
end
# TODO: 暂不考虑积雪模块

function Radiation(doy, hour, lat, Ta::T; kw...) where {T<:Real}
  rad = Radiation{T}(; kw...)
  update_rad!(rad, doy, hour, lat, Ta)
end

function update_rad!(rad::Radiation, doy, hour, lat, Ta::T=NaN) where {T<:Real}
  (; α_vis, α_nir, ϵ) = rad

  coszen = cal_coszen(doy, hour, lat)
  Rs_toa = cal_Rs_toa(doy, coszen) # Rs_toa
  Rs, Rs_dir, Rs_dif = partition_Rs(Rs_toa, coszen)

  # update Qa
  Rln_in = isnan(Ta) ? rad.Rln_in : _cal_Rli(Ta) # 重新更新Ta
  Qa = (1 - α_vis) * 0.5 * Rs + (1 - α_nir) * 0.5 * Rs + ϵ * Rln_in

  @pack! rad = Rs, Rs_dir, Rs_dif, Rs_toa, Rln_in, Qa, coszen
  rad
end

# bonan 2019, 案例
function _cal_Rli(Ta::Real)
  Θ = Ta + K0
  return (0.398e-05 * Θ^2.148) * σ * Θ^4 # 注意这里是气温，用气温推算的Rln_in
end

"""
# References
1. Gates, D.M. (1980), Clear sky atmospheric attenuation, Biophysical Ecology,
   110-115
"""
function partition_Rs(Rs_toa, coszen; τ_atm=0.5)
  oam = 1 / max(coszen, 0.04)
  Rs_dir = Rs_toa * τ_atm^oam # Clear sky direct beam
  Rs_dif = Rs_toa * (0.271 - 0.294 * τ_atm^oam) # Clear sky diffuse
  Rs = Rs_dif + Rs_dir # Clear sky total
  Rs, Rs_dir, Rs_dif
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


export _cal_Rli, update_rad!
