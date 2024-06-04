"""
    soil_thermal_properties(dz, mm_liq, mm_ice, Tsoil;
        soil_texture::Int=1, method="apparent-heat-capacity")

# Example
```julia
κ, cv = soil_thermal_properties(Tsoil, mm_liq, mm_ice, dz; soil_texture)
```

# Return
 - `κ`: Thermal conductivity, [W/m/K]
 - `cv`: Heat capacity of soil solids, [J/m3/K]
"""
function soil_thermal_properties!(κ, cv,
  Tsoil, SM_liq, SM_ice, dz;
  soil_texture::Int=1, method="apparent-heat-capacity")

  (; ρ_wat, ρ_ice, cv_wat, cv_ice, tk_wat, tk_ice) = physcon
  nsoil = length(dz)
  k = soil_texture
  TFRZ = K0
  cv_sol = 1.926 * 1e6 # Heat capacity of soil solids (J/m3/K)

  for i in 1:nsoil
    Θ_liq = SM_liq[i] / (ρ_wat * dz[i])  # [kg m-2] to [m3 m-3]
    Θ_ice = SM_ice[i] / (ρ_ice * dz[i])  # 
    fliq = Θ_liq / (Θ_liq + Θ_ice)       # Fraction of liq
    s = min((Θ_liq + Θ_ice) / Θ_S[k], 1) # Soil water relative to saturation

    # Dry thermal conductivity (W/m/K) from bulk density (kg/m3)
    bulk = 2700 * (1 - Θ_S[k])
    tk_dry = (0.135 * bulk + 64.7) / (2700 - 0.947 * bulk)

    # Thermal conductivity of quartz and soil solids (W/m/K)
    quartz = SAND[k] / 100     # Quartz fraction
    tko = quartz > 0.2 ? 2 : 3 # Thermal conductivity of other minerals (W/m/K)
    tk_quartz = 7.7
    tk_sol = tk_quartz^quartz * tko^(1 - quartz)

    # Saturated thermal conductivity (W/m/K) and unfrozen and frozen values
    tksat = tk_sol^(1 - Θ_S[k]) *
            tk_wat^(fliq * Θ_S[k]) *
            tk_ice^(Θ_S[k] - fliq * Θ_S[k])
    tksat_u = tk_sol^(1 - Θ_S[k]) * tk_wat^Θ_S[k]
    tksat_f = tk_sol^(1 - Θ_S[k]) * tk_ice^Θ_S[k]

    # Kersten number and unfrozen and frozen values
    if SAND[k] < 50
      ke_u = log10(max(s, 0.1)) + 1
    else
      ke_u = 0.7 * log10(max(s, 0.05)) + 1
    end
    ke_f = s
    ke = Tsoil[i] >= TFRZ ? ke_u : ke_f

    # Thermal conductivity (W/m/K) and unfrozen and frozen values
    κ[i] = (tksat - tk_dry) * ke + tk_dry
    κu = (tksat_u - tk_dry) * ke_u + tk_dry
    κf = (tksat_f - tk_dry) * ke_f + tk_dry

    # Heat capacity (J/m3/K) and unfrozen and frozen values
    cv[i] = (1 - Θ_S[k]) * cv_sol + cv_wat * Θ_liq + cv_ice * Θ_ice
    cvu = (1 - Θ_S[k]) * cv_sol + cv_wat * (Θ_liq + Θ_ice)
    cvf = (1 - Θ_S[k]) * cv_sol + cv_ice * (Θ_liq + Θ_ice)

    # Adjust heat capacity and thermal conductivity if using apparent heat capacity
    if method == "apparent-heat-capacity"
      tinc = 0.5 # Temperature range for freezing and thawing (K)

      if -tinc <= Tsoil[i] - TFRZ <= tinc
        # Heat of fusion (J/m3) - This is equivalent to ql = λ_fus * (mm_liq + mm_ice) / dz
        ql = λ_fus * (ρ_wat * Θ_liq + ρ_ice * Θ_ice)
        cv[i] = (cvf + cvu) / 2 + ql / (2 * tinc)
        κ[i] = κf + (κu - κf) * (Tsoil[i] - TFRZ + tinc) / (2 * tinc)
      elseif Tsoil[i] - TFRZ > tinc # unfrozen
        cv[i] = cvu
        κ[i] = κu
      elseif Tsoil[i] - TFRZ < -tinc # frozen
        cv[i] = cvf
        κ[i] = κf
      end
    end
  end
  κ, cv
end


export soil_thermal_properties!
