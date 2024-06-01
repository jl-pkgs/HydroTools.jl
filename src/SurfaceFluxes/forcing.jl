function RH2ea(RH::Real, Tair::Real=25.0)
  cal_es(Tair) * RH / 100
end

"""
    cal_ρ(Ta, ea, Pa) 
    
# Arguments
- `Ta`: Air temperature (degC)
- `ea`: Water vapor pressure (kPa)
- `Pa`: Air pressure (kPa)

# Return
- `ρₐ`: kg/m^3 (≈1.225 kg/m3)
"""
function cal_ρ(ea, Ta, Pa=atm)
  R = 8.31446         # Universal gas constant (J/K/mol)  
  ϵ = Mw / Md

  Ta += K0
  Rd = R / Md         # g mol-1
  ρₐ = (Pa - (1 - ϵ) * ea) / (Rd * Ta) # 抵消g的1e-3
  ρₐ
end

"""
    cal_θ(Pa_surf::FT, Ta_atm::FT) where {FT<:Real}

# Arguments
- `Pa_surf`: kPa
- `Ta_atm`: K
"""
function cal_θ(Pa_surf::FT, Ta_atm::FT) where {FT<:Real}
  m = Rd / (Cp * 1e6) # 0.283
  θ = Ta_atm * (Pa_surf / atm)^m
  θ
end

"""
    init_forcing(Ta, ea, Pa, z)

# Return
Meteorological variables at reference height z:
- `q`    : Specific humidity (kg/kg)
- `θ`    : Potential temperature (K)
- `θv`   : Virtual potential temperature (K)
- `ρₐ`   : Air density (kg/m3)
- `Mₐ`   : Molecular mass of air (kg/mol)
- `cpₐ`  : Specific heat of air at constant pressure, (J/mol/K)
- `ρ_mol`: Molar density (mol/m3)
"""
function init_forcing(Ta, ea, Pa, z)
  Ta += T0

  R = 8.31446         # Universal gas constant (J/K/mol)  
  ϵ = M_h2o / M_dry
  (; cp_w, cp_d) = physcon

  θ = Ta + Γd * z # potential temperature at reference height z
  q = ea2q(ea, Pa)

  θv = θ * (1 + 0.61 * q) # Tv

  ρ_mol = Pa / (R * Ta) # PV = nRT, p = ρ R T, 适用于任何气体
  ρₐ = ρ_mol * M_dry * (1 - (1 - ϵ) * ea / Pa) # kg m-3
  M_air = ρₐ / ρ_mol # kg mol-1
  
  cpₐ = cp_d * (1 + (cp_w / cp_d - 1) * q) # [J kg-1 K-1]
  @show cpₐ
  cpₐ = cpₐ * M_air # [J mol-1 K-1]
  return (;q, θ, θv, cpₐ, M_air, ρₐ, ρ_mol)
end


function cal_Rs_toa(lat, doy, hour)
  # Solar radiation at top of the atmosphere
  solcon = 1364 # W/m2

  decl = 23.45 * sin((284 + doy) / 365 * 2 * pi) * pi / 180
  hour_angle = 15 * (hour - 12) * pi / 180
  coszen = max(cos(lat) * cos(decl) * cos(hour_angle) + sin(lat) * sin(decl), 0)
  rv = 1 / sqrt(1 + 0.033 * cos(doy / 365 * 2 * pi))
  Rs_toa = solcon / rv^2 * coszen

  # Clear sky atmospheric attenuation: Gates, D.M. (1980) Biophysical Ecology, page 110, 115
  tau_atm = 0.5
  oam = 1 / max(coszen, 0.04)
  Rs_dir = Rs_toa * tau_atm^oam # Clear sky direct beam
  Rs_dif = Rs_toa * (0.271 - 0.294 * tau_atm^oam) # Clear sky diffuse
  Rs = Rs_dif + Rs_dir # Clear sky total

  return Rs_toa, Rs, Rs_dir, Rs_dif, coszen
end

export init_forcing, cal_Rs_toa
