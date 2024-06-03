export surface_fluxes

"""
Use Monin-Obukhov Similarity Theory (MOST) to obtain the Obukhov length (ζ).

This is the function to solve for the Obukhov length. For the current
estimate of the Obukhov length (x), calculate u*, T*, and q* and then
the new length (obu). The function value is the change in Obukhov length:
fx = x - obu.

# Examples
```julia
δζ = MOST(ζ, met, flux; z, d, z0m, z0c)
```
"""
function MOST(ζ0, met::Met, flux::Flux; z::Real, d::Real, z0m::Real, z0c::Real)
  k = 0.4             # von Karman constant
  g = 9.80665         # Gravitational acceleration (m/s2)
  # (; z0m, z0c, z, d) = param
  (; θ, e, u, Pa, M_air) = met
  (; θ_surf, e_surf) = flux
  ζ0 = abs(ζ0) <= 0.1 ? 0.1 : ζ0 # Prevent near-zero values
  z_minus_d = z - d

  # Evaluate psi for momentum at the reference height (zref-disp) and surface (z0m)
  Ψm = -Ψ_m_monin_obukhov(z0m / ζ0, z_minus_d / ζ0) # Eq 6.39
  Ψc = -Ψ_c_monin_obukhov(z0c / ζ0, z_minus_d / ζ0) # Eq 6.39

  zlog_m = log(z_minus_d / z0m)
  zlogc = log(z_minus_d / z0c)

  # Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)
  u₊ = u * k / (zlog_m + Ψm)                          # Eq. 6.40
  θ₊ = (θ - θ_surf) * k / (zlogc + Ψc)               # Eq. 6.43
  q₊ = (e - e_surf) / Pa * k / (zlogc + Ψc)          # Eq. 6.44
  θv₊ = θ₊ + 0.61 * θ * q₊ * (M_h2o / M_air)          # Eq. 6.33
  ζ = u₊^2 * θ / (k * g * θv₊)                        # Eq. 6.30

  @pack! flux = u₊, θ₊, q₊, ζ
  return ζ0 - ζ                                       # changes
end

function soil_temperature_delta(soil::Soil, df0::Real, f0::Real, snow_water::Real=0.0)
  (; dz, dt, κ, cv, Tsoil) = soil
  soil_temperature_delta(dz, dt, κ, cv, Tsoil, df0, f0, snow_water) # Tsoil_next, G_soil, G_snow
end

"""
# Arguments
- `fun`: most or rsl
  + `most`: Monin-Obukhov similarity theory
  + `rsl`: Roughness sublayer theory
"""
function surface_fluxes(met::Met, rad::Radiation, can::Canopy, soil::Soil, flux::Flux;
  param, snow_water=0.0, dt=1800.0)

  # Solve for the Obukhov length (m)
  function _most(ζ)
    (; z, d, z0m, z0c) = param
    MOST(ζ, met, flux; z, d, z0m, z0c) # δζ
  end
  ζ = root_hybrid(_most; tol=0.01, lb=100.0, ub=-100.0)
  @show ζ

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
  LW = ϵ * σ * Ts^4
  d_LW = 4 * ϵ * σ * Ts^3

  H = cpₐ * (θ_surf - θ) * g_ac # gH = g_ac
  d_H = cpₐ * g_ac

  es, d_es = satvap(θ_surf - K0) # ! bug may here
  LE = λ / Pa * (es - e) * gw * β_soil # miss a ϵ here
  d_LE = λ / Pa * d_es * gw * β_soil

  # Net energy flux into soil (W/m2) and temperature derivative (W/m2/K)
  f0 = Qa - LW - H - LE
  df0 = -d_LW - d_H - d_LE

  # solve Tsoil_next
  Tsoil_next, G_soil, G_snow = soil_temperature_delta(soil, df0, f0, snow_water)

  # Update surface fluxes for the change in surface temperature
  dtsrf = Tsoil_next[1] - Tsoil[1]
  LW += d_LW * dtsrf
  H += d_H * dtsrf
  LE += d_LE * dtsrf
  Rn = Qa - LW

  @pack! flux = g_ac, LE, H, G_soil, G_snow
  x = (; Qa, Rn, LE, H, G_soil, G_snow)
  @show x
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
