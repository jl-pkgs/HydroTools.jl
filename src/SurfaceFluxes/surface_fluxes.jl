"""
Calculate soil temperatures and surface fluxes. This routine uses the
current estimate of surface temperature and vapor pressure to solve
for the Obukhov length. This then provides the aerodynamic conductance,
which is used in the surface temperature and flux calculation.

Input
  dt                  ! Time step (s)

  forcvar.thref       ! Potential temperature at reference height (K)
  forcvar.uref        ! Wind speed at reference height (m/s)
  forcvar.eref        ! Vapor pressure at reference height (Pa)
  forcvar.pref        ! Atmospheric pressure (Pa)
  forcvar.cpair       ! Specific heat of air at constant pressure, at reference height (J/mol/K)
  ρ_mol      ! Molar density at reference height (mol/m3)

  surfvar.emiss       ! Surface emissivity
  surfvar.gcan        ! Canopy conductance (mol/m2/s)
  soilvar.method      ! Use excess heat or apparent heat capacity for phase change
  soilvar.nsoi        ! Number of soil layers
  soilvar.dz          ! Soil layer thickness (m)
  soilvar.cv          ! Heat capacity (J/m3/K)

  fluxvar.profiles    ! Use MOST or RSL for flux-profiles
  fluxvar.qa          ! Radiative forcing (W/m2)
  fluxvar.bucket      ! Use bucket model hydrology soil wetness factor

  bucket.snow_water   ! Snow water (kg H2O/m2)
  bucket.soil_water   ! Soil water (kg H2O/m2)
  bucket.soil_water_max ! Maximum soil water (kg H2O/m2)
  bucket.soil_beta_max  ! Soil water at which soil_beta = 1 (fraction of soil_water_max)

## Input/output
  Ta        ! Surface temperature (K)
  fluxvar.esrf        ! Surface vapor pressure (Pa)
  soilvar.tsoi        ! Soil temperature (K)
  soilvar.h2osoi_liq  ! Unfrozen water, liquid (kg H2O/m2)
  soilvar.h2osoi_ice  ! Frozen water, ice (kg H2O/m2)

## Output
  fluxvar.rnet        ! Net radiation (W/m2)
  fluxvar.lwrad       ! Emitted longwave radiation (W/m2)
  fluxvar.shflx       ! Sensible heat flux (W/m2)
  fluxvar.lhflx       ! Latent heat flux (W/m2)
  fluxvar.etflx       ! Evapotranspiration (mol H2O/m2/s)
  fluxvar.gsoi        ! Soil energy flux (W/m2)
  fluxvar.gsno        ! Snow melt energy flux (W/m2)
  fluxvar.gam         ! Aerodynamic conductance for momentum (mol/m2/s)
  fluxvar.gac         ! Aerodynamic conductance for scalars (mol/m2/s)
  soilvar.hfsoi       ! Soil phase change energy flux (W/m2)
  bucket.soil_beta    ! Soil wetness factor for evapotranspiration (-)
  bucket.snow_melt    ! Snow melt (kg H2O/m2/s)

## Output from Obukhov length calculation
  u₊       ! Friction velocity (m/s)
  t₊       ! Temperature scale (K)
  fluxvar.q₊       ! Water vapor scale (mol/mol)
  fluxvar.obu         ! Obukhov length (m)
  fluxvar.z0m         ! RSL only: Roughness length for momentum (m)
  fluxvar.z0c         ! RSL only: Roughness length for scalars (m)
  fluxvar.disp        ! RSL only: Displacement height (m)

Calculate the Obukhov length

Calculate the Obukhov length (obu) for the current surface temperature
and surface vapor pressure using Monin-Obukhov similarity theory or
Harman & Finnigan (2007, 2008) roughness sublayer (RSL) theory. Use the
functions "most" or "rsl" to iterate obu until the change in obu is less
than tol.
"""

"""
Use Monin-Obukhov similarity theory to obtain the Obukhov length (obu).

This is the function to solve for the Obukhov length. For the current
estimate of the Obukhov length (x), calculate u*, T*, and q* and then
the new length (obu). The function value is the change in Obukhov length:
fx = x - obu.
"""
function most(ζ, forcvar, fluxvar, param)
  k = 0.4             # von Karman constant
  g = 9.80665         # Gravitational acceleration (m/s2)

  (; z0m, z0c, z, d) = param
  (; θ_ref, e_ref, u_ref, Pa_ref, M_air) = forcvar
  (; θ_surf, e_surf) = fluxvar
  # Prevent near-zero values of the Obukhov length
  ζ = abs(ζ) <= 0.1 ? 0.1 : ζ

  # Calculate z-d at the reference height, because this is used many times
  z_minus_d = z - d

  # Evaluate psi for momentum at the reference height (zref-disp) and surface (z0m)
  Ψm = -Ψ_m_monin_obukhov(z0m / ζ, z_minus_d / ζ) # Eq 6.39
  Ψc = -Ψ_c_monin_obukhov(z0c / ζ, z_minus_d / ζ) # Eq 6.39

  zlog_m = log(z_minus_d / z0m)
  zlog_c = log(z_minus_d / z0c)

  # Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)
  u₊ = u_ref * k / (zlog_m + Ψm)                          # Eq. 6.40
  θ₊ = (θ_ref - θ_surf) * k / (zlog_c + Ψc)               # Eq. 6.43
  q₊ = (e_ref - e_surf) / Pa_ref * k / (zlog_c + Ψc)          # Eq. 6.44
  θv₊ = θ₊ + 0.61 * θ_ref * q₊ * (M_h2o / M_air)          # Eq. 6.33

  ζ_next = u₊^2 * θ_ref / (k * g * θv₊)  # Obukhov length (m), Eq 6.30
  return ζ - ζ_next                       # changes
end


"""
# Arguments
- `fun`: most or rsl
  + `most`: Monin-Obukhov similarity theory
  + `rsl`: Roughness sublayer theory
"""
function surface_fluxes(physcon, forcvar, surfvar, soilvar, fluxvar, bucket, dt)
  obu_0 = 100  # Initial estimate for Obukhov length (m)
  obu_1 = -100  # Initial estimate for Obukhov length (m)
  tol = 0.01  # Accuracy tolerance for Obukhov length (m)

  # Solve for the Obukhov length
  fluxvar, _ = root_hybrid(most, obu_0, obu_1, tol, forcvar, surfvar, fluxvar)
  # ? solve for the Obukhov length

  (; λ_sub, λ_vap, M_h2o) = physcon
  (; snow_water, soil_water, soil_beta_max, soil_water_max) = bucket
  (; gcan) = surfvar
  (; gac, u₊, t₊) = fluxvar
  (; θ_ref, e_ref, Pa_ref) = forcvar

  # Aerodynamic conductances for momentum (gam) and scalars (gac) (mol/m2/s)
  # gam = ρ_mol * u₊^2 / uref # not used
  gac = ρ_mol * u₊ * t₊ / (θ_ref - θ_surf) # Eq 6.7 and 6.15
  gw = 1 / (1 / gcan + 1 / gac) # Surface conductance for water vapor (mol/m2/s)

  # Latent heat of vaporization or sublimation (J/mol)
  λ = snow_water > 0 ? λ_sub : λ_vap
  λ *= M_h2o # [J kg-1 * kg mol-1] = J mol-1

  # Soil wetness factor for evapotranspiration
  β_soil = bucket == "no_bucket" ? 1 : min(soil_water / (soil_beta_max * soil_water_max), 1)

  Ta = θ_surf
  # Emitted longwave radiation (W/m2) and temperature derivative (W/m2/K)
  LW = ϵ * σ * Ta^4
  d_LW = 4 * ϵ * σ * Ta^3
  
  H = cp_air * (θ_surf - θ_ref) * gH # gH = gac
  d_H = cp_air * gH

  es, d_es = satvap(θ_surf - K0) # ! bug may here
  LE = λ / Pa_ref * (es - e_ref) * gw * β_soil # miss a ϵ here
  d_LE = λ / Pa_ref * d_es * gw * β_soil

  # Net energy flux into soil (W/m2) and temperature derivative (W/m2/K)
  f0 = Qa - LW - H - LE
  df0 = -d_LW - d_H - d_LE

  # Update soil temperatures
  Tsoil_next, G_soil, G_snow = soil_temperature_delta(dz, dt, κ, cv, Tsoil_cur,
    df0, f0, snow_water)

  # Update surface fluxes for the change in surface temperature
  dtsrf = Tsoil_next[1] - Tsoil[1]
  LW += d_LW * dtsrf
  H += d_H * dtsrf
  LE += d_LE * dtsrf
  Rn = Qa - LW

  # Error check
  err = Rn - LE - H - G_soil - G_snow
  if abs(err) > 1e-06
    println("err = ", err)
    x = (; Qa, Rn, LE, H, G_soil, G_snow)
    @show x
    error("surface temperature error")
  end

  ET = LE / λ # Evapotranspiration (mol H2O/m2/s), [W m-2 / (J mol-1)]
  e_surf = (e_ref / Pa_ref + ET / gac) * Pa_ref
  # Surface vapor pressure is diagnosed from evaporative flux
  
  # Phase change for soil layers undergoing freezing of thawing
  # if soilvar.method == "apparent-heat-capacity"
  hfsoi = 0
  # elseif soilvar.method == "excess-heat"
  # phase_change!(physcon, soilvar, dt)
  # end

  # Check for energy conservation
  edif = 0
  for i in 1:nsoil
    edif += cv[i] * dz[i] * (Tsoil_next[i] - Tsoil[i]) / dt
  end
  err = edif - G_soil - hfsoi
  abs(err) > 1e-03 && error("Soil temperature energy conservation error")

  # Surface temperature is the first soil layer
  Ts = Tsoil_next[1]
end

# ea2ρ(ea) = ϵ * ea / Pa * ρₐ
