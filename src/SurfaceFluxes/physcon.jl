module physcon

# vkc    = 0.4             # von Karman constant
# grav   = 9.80665         # Gravitational acceleration (m/s2)
tfrz   = 273.15          # Freezing point of water (K)
σ      = 5.67e-08        # Stefan-Boltzmann constant (W/m2/K4)

M_dry  = 28.97 / 1000    # Molecular mass of dry air (kg/mol)
M_h2o  = 18.02 / 1000    # Molecular mass of water (kg/mol)

rgas   = 8.31446         # Universal gas constant (J/K/mol)
ρ_wat  = 1000.0          # Density of water (kg/m3)
ρ_ice  = 917.0           # Density of ice (kg/m3)
cp_wat = 4188.0          # Specific heat of water (J/kg/K)
cp_ice = 2117.27         # Specific heat ice (J/kg/K)
cv_wat = cp_wat * ρ_wat  # Heat capacity of water (J/m3/K)
cv_ice = cp_ice * ρ_ice  # Heat capacity of ice (J/m3/K)
tk_wat = 0.57            # Thermal conductivity of water (W/m/K)
tk_ice = 2.29            # Thermal conductivity of ice (W/m/K)
λ_fus  = 0.3337 * 1e6    # Heat of fusion for water at 0 C (J/kg)
λ_vap  = 2.501 * 1e6     # Latent heat of evaporation (J/kg)
λ_sub  = λ_fus + λ_vap   # Latent heat of sublimation (J/kg)
solcon = 1367            # W/m2

end
