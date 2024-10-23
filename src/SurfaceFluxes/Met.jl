export Met, update_met!

"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Met{T<:AbstractFloat}
  "reference height (m)"
  z::T = NaN # m
  "speed (m/s)"
  u::T = NaN
  "air temperature (C)"
  Ta::T = 20.0 # degC
  "specific humidity (kg/kg)"
  q::T = NaN
  "pressure (Pa)"
  Pa::T = 101325.0
  "potential temperature (C)"
  θ::T = Ta + Γd * z
  "virtual potential temperature (C)"
  θv::T = θ * (1 + 0.61 * q)
  "vapor pressure (Pa)"
  e::T = q2ea(q, Pa)
  "rainfall (kg/m2, mm)"
  rain::T = NaN
  "snowfall (kg/m2, mm)"
  snow::T = NaN
  "Specific heat of air at constant pressure, (J/mol/K)"
  cpₐ::T = NaN
  "Molecular mass of air (kg/mol)"
  M_air::T = NaN
  "Air density (kg/m3)"
  ρₐ::T = NaN
  "Molar density (mol/m3)"
  ρ_mol::T = NaN
end

"""
    update_met!(met, Ta, e, u, rain, snow, Pa)

# Notes
提前设定`z`，Reference height (m)，用于计算位温

# Arguments
- `e`  : Vapor pressure (Pa)
- `Pa` : Air pressure (Pa)，用于计算q, ρ_mol, ρₐ，影响应该不大

# Return
Meteorological variables at reference height z:
- `q`    : Specific humidity (kg/kg)
- `θ`    : Potential temperature (K)
- `θv`   : Virtual potential temperature (K)
- `cpₐ`  : Specific heat of air at constant pressure, (J/mol/K)
- `M_air`: Molecular mass of air (kg/mol)
- `ρₐ`   : Air density (kg/m3)
- `ρ_mol`: Molar density (mol/m3)
"""
function update_met!(met::Met, Ta::T, e::T, u::T, rain::T=0, snow::T=0, Pa::T=atm * 1000) where {T<:Real}
  (; z) = met
  Ta += K0
  # https://github.com/jl-pkgs/bonanmodeling/blob/master/sp_07_Surface%20Energy%20Fluxes/physcon.m
  ϵ = M_h2o / M_dry
  R = 8.31446        # Universal gas constant (J/K/mol)  
  cp_d = 1005.0      # Specific heat of `dry air` at constant pressure (J/kg/K)
  cp_w = 1846.0      # Specific heat of `water vapor` at constant pressure (J/kg/K)

  θ = Ta + Γd * z # potential temperature at reference height z
  q = ea2q(e, Pa)

  θv = θ * (1 + 0.61 * q) # Tv

  ρ_mol = Pa / (R * Ta) # PV = nRT, p = ρ R T, 适用于任何气体, mol m-3
  ρₐ = ρ_mol * M_dry * (1 - (1 - ϵ) * e / Pa) # kg m-3
  M_air = ρₐ / ρ_mol # kg mol-1

  cpₐ = cp_d * (1 + (cp_w / cp_d - 1) * q) # [J kg-1 K-1]
  cpₐ = cpₐ * M_air # [J mol-1 K-1]

  @pack! met = q, θ, θv, cpₐ, M_air, ρₐ, ρ_mol,
    Pa, e, u, rain, snow
  met
end
