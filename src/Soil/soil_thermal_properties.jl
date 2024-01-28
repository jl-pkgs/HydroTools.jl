"""
# Arguments

- `method`: method of phase change

"""
function soil_thermal_properties(physcon, soilvar; method="excess-heat")

  ρ_wat = 1000.0                       # Density of water (kg/m3)
  ρ_ice = 917.0                        # Density of ice (kg/m3)

  _c_wat = 4188.0                         # Specific heat of water (J/kg/K)
  _c_ice = 2117.27                        # Specific heat of ice (J/kg/K)

  cv_wat = _c_wat * ρ_wat # Heat capacity of water (J/m3/K)
  cv_ice = _c_ice * ρ_ice # Heat capacity of ice (J/m3/K)

  tfrz = 273.15                         # Freezing point of water [k]

  ## --- Physical constants in physcon structure
  κ_wat = 0.57                          # Thermal conductivity of water (W/m/K)
  κ_ice = 2.29                          # Thermal conductivity of ice (W/m/K)
  physcon.hfus = 0.3337e6               # Heat of fusion for water at 0 C (J/kg)

  ## --- Model run control parameters
  soil_texture = 1          # Soil texture class: sand

  @unpack dz, h2osoi_liq, h2osoi_ice, tsoi = soilvar

  for i = 1:soilvar.nsoi
    k = soil_texture

    # --- Volumetric soil water and ice
    Θ_liq = h2osoi_liq[i] / (ρ_wat * dz[i])
    Θ_ice = h2osoi_ice[i] / (ρ_ice * dz[i])

    # Fraction of total volume that is liquid water
    fᵤ = Θ_liq / (Θ_liq + Θ_ice)

    # Soil water relative to saturation
    s = min((Θ_liq + Θ_ice) / Θ_S[k], 1)

    # --- Dry thermal conductivity (W/m/K) from 
    ρ_b = 2700 * (1 - Θ_S[k]) # density of soil soliads, bulk density (kg/m3)
    κ_dry = (0.135 * ρ_b + 64.7) / (2700 - 0.947 * ρ_b)  # Eq. 5.27

    ## --- Soil solids thermal conducitivty (W/m/K)
    q = SAND[k] / 100        # Quartz fraction
    κ_o = q > 0.2 ? 2.0 : 3.0 # Thermal conductivity of other minerals (W/m/K)
    κ_q = 7.7                # Thermal conductivity of q (W/m/K)

    # Thermal conductivity of soil solids (W/m/K)
    κ_sol = κ_q^q * κ_o^(1 - q)  # Eq. 5.31

    # --- Saturated thermal conductivity (W/m/K) and unfrozen and frozen values
    κ_sat = κ_sol^(1 - Θ_S[k]) * κ_wat^(fᵤ * Θ_S[k]) * κ_ice^((1 - fᵤ) * Θ_S[k]) # Eq. 5.30

    κ_sat_u = κ_sol^(1 - Θ_S[k]) * κ_wat^Θ_S[k] # Eq. 5.28
    κ_sat_f = κ_sol^(1 - Θ_S[k]) * κ_ice^Θ_S[k] # Eq. 5.29

    # --- Kersten number and unfrozen and frozen values
    if (sand[k] < 50)
      ke_u = log10(max(s, 0.1)) + 1
    else
      ke_u = 0.7 * log10(max(s, 0.05)) + 1
    end
    ke_f = s

    ke = tsoi[i] >= tfrz ? ke_u : ke_f

    # --- Thermal conductivity (W/m/K) and unfrozen and frozen values
    soilvar.κ[i] = (κ_sat - κ_dry) * ke + κ_dry
    κ_u = (κ_sat_u - κ_dry) * ke_u + κ_dry
    κ_f = (κ_sat_f - κ_dry) * ke_f + κ_dry

    # --- Heat capacity of soil solids (J/m3/K)
    cv_sol = 1.926e6

    # --- Heat capacity (J/m3/K) and unfrozen and frozen values
    soilvar.cv[i] = (1 - Θ_S[k]) * cv_sol + cv_wat * Θ_liq + cv_ice * Θ_ice
    cvu = (1 - Θ_S[k]) * cv_sol + cv_wat * (Θ_liq + Θ_ice)
    cvf = (1 - Θ_S[k]) * cv_sol + cv_ice * (Θ_liq + Θ_ice)

    # --- Adjust heat capacity and thermal conductivity if using apparent heat capacity
    #  switch soilvar.method
    #     case "apparent-heat-capacity"
    #     # Temperature range for freezing and thawing [k]
    #     tinc = 0.5;

    #     # Heat of fusion (J/m3) - This is equivalent to ql = hfus * (h2osoi_liq + h2osoi_ice) / dz
    #     ql = physcon.hfus * (ρ_wat * Θ_liq + ρ_ice * Θ_ice);

    #     # Heat capacity and thermal conductivity
    #     if (soilvar.tsoi[i] > tfrz+tinc)
    #        soilvar.cv[i] = cvu;
    #        soilvar.κ_[i] = κ_u;
    #     end

    #     if (soilvar.tsoi[i] >= tfrz-tinc & soilvar.tsoi[i] <= tfrz+tinc)
    #        soilvar.cv[i] = (cvf + cvu) / 2 + ql / (2 * tinc);
    #        soilvar.κ_[i] = κ_f + (κ_u - κ_f) * (soilvar.tsoi[i] - tfrz + tinc) / (2 * tinc);
    #     end

    #     if (soilvar.tsoi[i] < tfrz-tinc)
    #        soilvar.cv[i] = cvf;
    #        soilvar.κ_[i] = κ_f;
    #     end
    #  end
  end
end
