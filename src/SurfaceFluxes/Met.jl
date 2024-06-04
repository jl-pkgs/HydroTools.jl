"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Met{T<:AbstractFloat}
  "reference height (m)"
  z::T = 0.0 # m
  "speed (m/s)"
  u::T = T(NaN)
  "air temperature (C)"
  Ta::T = 20.0 # degC
  "specific humidity (kg/kg)"
  q::T = 0.0
  "pressure (Pa)"
  Pa::T = 101325.0
  "potential temperature (C)"
  θ::T = Ta + Γd * z
  "virtual potential temperature (C)"
  θv::T = θ * (1 + 0.61 * q)
  "vapor pressure (Pa)"
  e::T = q2ea(q, Pa)
  "rainfall (kg/m2, mm)"
  rain::T = T(NaN)
  "snowfall (kg/m2, mm)"
  snow::T = T(NaN)
  "Specific heat of air at constant pressure, (J/mol/K)"
  cpₐ::T = T(NaN)
  "Molecular mass of air (kg/mol)"
  M_air::T = T(NaN)
  "Air density (kg/m3)"
  ρₐ::T = T(NaN)
  "Molar density (mol/m3)"
  ρ_mol::T = T(NaN)
end

"""
    init_forcing(Ta, ea, Pa, z)

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
function init_met(Ta, ea, Pa, z)
  Ta += K0

  # https://github.com/jl-pkgs/bonanmodeling/blob/master/sp_07_Surface%20Energy%20Fluxes/physcon.m
  ϵ = M_h2o / M_dry
  R = 8.31446        # Universal gas constant (J/K/mol)  
  cp_d = 1005.0      # Specific heat of `dry air` at constant pressure (J/kg/K)
  cp_w = 1846.0      # Specific heat of `water vapor` at constant pressure (J/kg/K)

  θ = Ta + Γd * z # potential temperature at reference height z
  q = ea2q(ea, Pa)

  θv = θ * (1 + 0.61 * q) # Tv

  ρ_mol = Pa / (R * Ta) # PV = nRT, p = ρ R T, 适用于任何气体, mol m-3
  ρₐ = ρ_mol * M_dry * (1 - (1 - ϵ) * ea / Pa) # kg m-3
  M_air = ρₐ / ρ_mol # kg mol-1

  cpₐ = cp_d * (1 + (cp_w / cp_d - 1) * q) # [J kg-1 K-1]
  cpₐ = cpₐ * M_air # [J mol-1 K-1]
  return (; q, θ, θv, cpₐ, M_air, ρₐ, ρ_mol)
end

function Met(Ta, ea, Pa, z; kw...)
  r = init_met(Ta, ea, Pa, z)
  Met{Float64}(; r..., z, kw...)
end
