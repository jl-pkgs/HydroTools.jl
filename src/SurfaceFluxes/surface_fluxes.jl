"""
# Arguments
- `fun`: most or rsl
  + `most`: Monin-Obukhov similarity theory
  + `rsl`: Roughness sublayer theory
"""
function surface_fluxes(met::Met, rad::Radiation, can::Canopy, soil::Soil, flux::Flux;
  param, snow_water=0.0)
  
  # Solve for the Obukhov length (m)
  most(x) = MOST(x, met, flux; param...)
  ζ = root_hybrid(most; tol=0.01, lb=100.0, ub=-100.0) # fill! flux; 有时这里会求解失败

  # (; snow_water, soil_water, soil_beta_max, soil_water_max) = bucket
  (; θ, e, Pa, ρ_mol, cpₐ) = met
  (; ϵ, Qa) = rad
  (; gc) = can
  (; θ_surf, e_surf, g_ac, u₊, θ₊) = flux
  (; Tsoil) = soil

  # Aerodynamic conductances for momentum (gam) and scalars (g_ac) (mol/m2/s)
  # gam = ρ_mol * u₊^2 / uref # not used
  g_ac = ρ_mol * u₊ * θ₊ / (θ - θ_surf) # Eq 6.7 and 6.15
  gw = 1 / (1 / gc + 1 / g_ac) # Surface conductance for water vapor (mol/m2/s)

  # Latent heat of vaporization or sublimation (J/mol)
  λ = snow_water > 0 ? λ_sub : λ_vap
  λ *= M_h2o # [J kg-1 * kg mol-1] = J mol-1

  # Soil wetness factor for evapotranspiration
  # β_soil = bucket == "no_bucket" ? 1 : min(soil_water / (soil_beta_max * soil_water_max), 1)
  β_soil = 1

  Ts = θ_surf
  # Emitted longwave radiation (W/m2) and temperature derivative (W/m2/K)
  LWout = ϵ * σ * Ts^4 # 向外的部分
  d_LWout = 4 * ϵ * σ * Ts^3

  H = cpₐ * (θ_surf - θ) * g_ac # gH = g_ac
  d_H = cpₐ * g_ac

  es, d_es = satvap(θ_surf - K0) # ! bug may here
  LE = λ / Pa * (es - e) * gw * β_soil # miss a ϵ here
  d_LE = λ / Pa * d_es * gw * β_soil

  # Net energy flux into soil (W/m2) and temperature derivative (W/m2/K)
  f0 = Qa - LWout - H - LE
  df0 = -d_LWout - d_H - d_LE

  # solve Tsoil_next
  Tsoil_next, G_soil, G_snow = soil_temperature_delta(soil, df0, f0, snow_water)

  # Update surface fluxes for the change in surface temperature
  dtsrf = Tsoil_next[1] - Tsoil[1]
  LWout += d_LWout * dtsrf
  H += d_H * dtsrf
  LE += d_LE * dtsrf
  Rn = Qa - LWout

  @pack! flux = g_ac, Qa, LWout, Rn, LE, H, G_soil, G_snow
  flux
end

# # ea2ρ(ea) = ϵ * ea / Pa * ρₐ

# # Error check
# err = Rn - LE - H - G_soil - G_snow
# if abs(err) > 1e-06
#   println("err = ", err)
#   x = (; Qa, Rn, LE, H, G_soil, G_snow)
#   @show x
#   error("surface temperature error")
# end

# ET = LE / λ # Evapotranspiration, [W m-2 / (J mol-1)] to [mol H2O m-2 s-1]
# flux.e_surf = (e / Pa + ET / g_ac) * Pa
# # Surface vapor pressure is diagnosed from evaporative flux

# # Phase change for soil layers undergoing freezing of thawing
# # if soilvar.method == "apparent-heat-capacity"
# hfsoi = 0
# # elseif soilvar.method == "excess-heat"
# # phase_change!(physcon, soilvar, dt)
# # end

# # Check for energy conservation
# edif = 0
# for i in 1:nsoil
#   edif += cv[i] * dz[i] * (Tsoil_next[i] - Tsoil[i]) / dt
# end
# err = edif - G_soil - hfsoi
# abs(err) > 1e-03 && error("Soil temperature energy conservation error")

# # Surface temperature is the first soil layer
# # Ts = Tsoil_next[1]
