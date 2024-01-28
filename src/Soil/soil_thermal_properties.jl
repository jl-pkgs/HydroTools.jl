"""
    soil_thermal_properties(dz::AbstractVector, Tsoil::AbstractVector,
        m_liq::AbstractVector, m_ice::AbstractVector;
        soil_texture::Integer=1, method="excess-heat")

# Arguments

- `dz`: the thickness of each soil layer (m)
- `m_liq`: Unfrozen water, liquid (kg H2O/m2)
- `m_ice`: Frozen water, ice (kg H2O/m2)
- `Tsoil`: Soil temperature of each soil layer (K)

- `method`: method of phase change
- `soil_texture`: 
  + `1`: sand

# Return

- `κ` : thermal conductivity, W/m/K
- `cv`: heat capacity, J/m3/K
"""
function soil_thermal_properties(dz::AbstractVector, Tsoil::AbstractVector,
  m_liq::AbstractVector, m_ice::AbstractVector;
  soil_texture::Integer=1, 
  method="excess-heat")
  
  _c_wat = 4188.0                         # Specific heat of water (J/kg/K)
  _c_ice = 2117.27                        # Specific heat of ice (J/kg/K)

  cv_wat = _c_wat * ρ_wat # Heat capacity of water (J/m3/K)
  cv_ice = _c_ice * ρ_ice # Heat capacity of ice (J/m3/K)

  tfrz = 273.15                         # Freezing point of water [k]

  ## --- Physical constants in physcon structure
  κ_wat = 0.57                          # Thermal conductivity of water (W/m/K)
  κ_ice = 2.29                          # Thermal conductivity of ice (W/m/K)
  λ_fus = 0.3337e6                      # Heat of fusion for water at 0 C (J/kg)

  ## --- Model run control parameters
  n = length(dz)
  κ = zeros(n)
  cv = zeros(n)

  k = soil_texture
  for i = 1:n
    # --- Volumetric soil water and ice
    Θ_liq = m_liq[i] / (ρ_wat * dz[i])
    Θ_ice = m_ice[i] / (ρ_ice * dz[i])

    # Fraction of total volume that is liquid water
    fᵤ = Θ_liq / (Θ_liq + Θ_ice)

    # --- Dry thermal conductivity (W/m/K) from 
    ρ_b = 2700 * (1 - Θ_S[k]) # density of soil soliads, bulk density (kg/m3)
    κ_dry = (0.135 * ρ_b + 64.7) / (2700 - 0.947 * ρ_b)  # Eq. 5.27

    # --- Kersten number and unfrozen and frozen values
    S_e = min((Θ_liq + Θ_ice) / Θ_S[k], 1) # Soil water relative to saturation
    Ke_f = S_e

    if (SAND[k] < 50)
      Ke_u = 1 + log10(max(S_e, 0.1))
    else
      Ke_u = 1 + 0.7 * log10(max(S_e, 0.05))
    end
    Ke = Tsoil[i] >= tfrz ? Ke_u : Ke_f

    ## --- Soil solids thermal conducitivty (W/m/K)
    q = SAND[k] / 100         # Quartz fraction
    κ_o = q > 0.2 ? 2.0 : 3.0 # Thermal conductivity of other minerals (W/m/K)
    κ_q = 7.7                 # Thermal conductivity of q (W/m/K)

    # Thermal conductivity of soil solids (W/m/K)
    κ_sol = κ_q^q * κ_o^(1 - q)  # Eq. 5.31

    # --- Saturated thermal conductivity (W/m/K) and unfrozen and frozen values
    κ_sat = κ_sol^(1 - Θ_S[k]) * κ_wat^(fᵤ * Θ_S[k]) * κ_ice^((1 - fᵤ) * Θ_S[k]) # Eq. 5.30

    κ_sat_u = κ_sol^(1 - Θ_S[k]) * κ_wat^Θ_S[k] # Eq. 5.28
    κ_sat_f = κ_sol^(1 - Θ_S[k]) * κ_ice^Θ_S[k] # Eq. 5.29

    # --- Thermal conductivity (W/m/K) and unfrozen and frozen values
    κ[i] = (κ_sat - κ_dry) * Ke + κ_dry

    κ_u = (κ_sat_u - κ_dry) * Ke_u + κ_dry
    κ_f = (κ_sat_f - κ_dry) * Ke_f + κ_dry

    ## --- Heat capacity of soil solids (J/m3/K)
    cv_sol = 1.926e6

    # --- Heat capacity (J/m3/K) and unfrozen and frozen values
    cv[i] = (1 - Θ_S[k]) * cv_sol + cv_wat * Θ_liq + cv_ice * Θ_ice # Eq. 5.32

    cv_u = (1 - Θ_S[k]) * cv_sol + cv_wat * (Θ_liq + Θ_ice)
    cv_f = (1 - Θ_S[k]) * cv_sol + cv_ice * (Θ_liq + Θ_ice)

    # --- Adjust heat capacity and thermal conductivity if using apparent heat capacity
    if method == "apparent-heat-capacity"
      # 这里考虑了结冰和融化的过程
      tinc = 0.5 # Temperature range for freezing and thawing [k]
      # Heat of fusion (J/m3), equivalent to 
      # ql = λ_fus * (m_liq + m_ice) / dz
      ql = λ_fus * (ρ_wat * Θ_liq + ρ_ice * Θ_ice)

      # Heat capacity and thermal conductivity, Eq. 5.39
      if Tsoil[i] > tfrz + tinc
        cv[i] = cv_u
        κ[i] = κ_u
      elseif tfrz - tinc <= Tsoil[i] <= tfrz + tinc
        cv[i] = (cv_f + cv_u) / 2 + ql / (2 * tinc)
        κ[i] = κ_f + (κ_u - κ_f) * (Tsoil[i] - tfrz + tinc) / (2 * tinc)
      elseif Tsoil[i] < tfrz - tinc
        cv[i] = cv_f
        κ[i] = κ_f
      end
    end
  end
  κ, cv
end
