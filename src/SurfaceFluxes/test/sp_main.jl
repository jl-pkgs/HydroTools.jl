include("physcon.jl")
include("satvap.jl")

j = 1:nday
i = 1:ntime

i, j = 20, 1
hour = i * (dt / 86400 * 24)

Ta = tmean + 0.5 * trange * sin(2 * pi / 24 * (hour - 8)) + K0
es, d_es = satvap(Ta - K0)
ea = es * RH / 100

cpₐ, T_pot, Tv, mm_air, ρ_air, ρ_mol = cal_Cp(Ta, ea, Pa, z)

Rs_toa, Rs, Rs_dir, Rs_dif, coszen = cal_Rs_toa(lat, doy, hour)

σ = physcon.sigma
Rln_in = (0.398e-05 * Ta^2.148) * σ * Ta^4;

# Maximum soil water (kg H2O/m2)
bucket = Dict(
  "soil_water_max" => 150.0, 
  "soil_beta_max" => 0.75, 
  "snow_water" => 0, 
  "soil_water" => 150.0) # soil_water_max

fsno = bucket["snow_water"] / (bucket["snow_water"] + snow_mask)
alb_eff_vis = α_surf.vis * (1 - fsno) + α_snow.vis * fsno
alb_eff_nir = α_surf.nir * (1 - fsno) + α_snow.nir * fsno

# Radiative forcing: absorbed solar + incident longwave. This partitions
# solar radiation into 50% visible and 50% near-infrared wavebands.
Qa = (1 - alb_eff_vis) * 0.5 * Rs +
  (1 - alb_eff_nir) * 0.5 * Rs + ϵ * Rln_in

# Canopy conductance (mol/m2/s) - use a weighted average of sunlit and shaded leaves
gcan_min = 0.05  # Minimum conductance (mol/m2/s)
gcan_max = 0.2   # Maximum conductance (mol/m2/s)

ext = 0.5 / max(coszen, 0.0001)  # Light extinction coefficient
fsun = (1 - exp(-ext * LAI)) / (ext * LAI)  # Sunlit fraction of canopy

if (Rs > 0)
  gcan = (fsun * gcan_max + (1 - fsun) * gcan_min) * LAI
else
  gcan = gcan_min * LAI
end

# % Thermal conductivity and heat capacity
κ, cv = soil_thermal_properties(Tsoil, mm_liq, mm_ice, dz; soil_texture)

# Calculate the soil temperatures and surface fluxes
fluxvar, soilvar, bucket = surface_fluxes(forcvar, surfvar, soilvar, fluxvar, bucket, dt)

# Rainfall to equal evaporative loss (kg H2O/m2/s)
# forcvar.rain = fluxvar.etflx * physcon.mmh2o
forcvar.rain = 0

# Bucket model hydrology
bucket = bucket_hydrology(physcon, forcvar, fluxvar, bucket, dt)

# Save data for graphics
xhour[i] = hour
ytsrf[i] = fluxvar.tsrf - physcon.tfrz
ytref[i] = forcvar.tref - physcon.tfrz
yrnet[i] = fluxvar.rnet
yshflx[i] = fluxvar.shflx
ylhflx[i] = fluxvar.lhflx
ygsoi[i] = fluxvar.gsoi
yustar[i] = fluxvar.ustar
ygac[i] = fluxvar.gac * 100
ygcan[i] = surfvar.gcan * 100
