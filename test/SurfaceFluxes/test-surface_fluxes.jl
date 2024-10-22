using HydroTools, Test
using HydroTools.SurfaceFluxes

begin
  z = 30.0 # Reference height (m)
  hc = 20.0 # canopy height
  d = 0.67 * hc # Zero-plane displacement height (m)
  z0m = 0.13 * hc # Roughness length for momentum (m)
  z0c = 0.10 * z0m # Roughness length for scalars (m)
  param = (; z, d, z0m, z0c)
end

begin
  day = 182
  j = 20
  hour = j * dt / 86400 * 24
  lat = 40.0 * pi / 180  # Latitude (degrees -> radians) for solar radiation
  dt = 1800 # second

  ## Meteorological forcing 
  # Ta = 20.0
  Tmean = 25.0
  Ta = gen_Ta(hour)  # Air temperature (C)
  RH = 70.0  # Relative humidity (%)
  Pa = atm * 1e3
  es, d_es = satvap(Ta) # Pa
  ea = es * RH / 100
  z = 30.0
  met = Met(Ta, ea, Pa, z; rain=0, snow=0, u=3.0)

  ## Flux
  Ts = Tmean + K0
  es, d_es = satvap(Ts - K0)
  flux = Flux{Float64}(; θ_surf=Ts, e_surf=es)

  ## Radiation and Canopy
  Θ = Ta + K0
  Rln_in = (0.398e-05 * Θ^2.148) * σ * Θ^4
  rad = Radiation(day, hour, lat; Rln_in)

  coszen = cal_coszen(day, hour, lat)
  can = Canopy{Float64}(; LAI=5.0, coszen)

  ## Soil
  dz = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058]
  soil = Soil(dz)
  init_soil!(soil)
  soil
  surface_fluxes!(flux, met, rad, can, soil; param)
end

@testset "surface_fluxes" begin
  @test flux.Rn == flux.Qa - flux.LWout
  @test flux.Rn ≈ flux.H + flux.LE + flux.G_soil + flux.G_snow
  @test flux.H ≈ 44.301017430068185
  @test flux.LE ≈ 218.35986837001974
end

# Rainfall to equal evaporative loss (kg H2O/m2/s)
# forcvar.rain = fluxvar.etflx * physcon.mmh2o
# forcvar.rain = 0
# Bucket model hydrology
# bucket = bucket_hydrology(physcon, forcvar, fluxvar, bucket, dt)
# Ts = soil.Tsoil[1]
# @with_kw struct Flux{T<:AbstractFloat}
#   t_surf::T = Ts
#   e_surf::T = cal_es(Ts)
# end
