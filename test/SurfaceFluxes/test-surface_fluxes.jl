using HydroTools, Test
using HydroTools.SurfaceFluxes

function cal_Rln_in(Ta::Real)
  θ = Ta + K0
  return (0.398e-05 * θ^2.148) * σ * θ^4 # 注意这里是气温，用气温推算的Rln_in
end


begin
  z = 30.0 # Reference height (m)
  hc = 20.0 # canopy height
  d = 0.67 * hc # Zero-plane displacement height (m)
  z0m = 0.13 * hc # Roughness length for momentum (m)
  z0c = 0.10 * z0m # Roughness length for scalars (m)
  param = (; z, d, z0m, z0c)

  # 气象条件
  rain = 0.0
  snow = 0.0
  u = 3.0
  Pa = atm * 1e3
  RH = 70.0  # Relative humidity (%)
  Tmean = 25.0
end

met = Met(; z=30.0) # 首先设定高度

begin
  dt = 1800 # second
  day = 182
  j = 20
  hour = j * dt / 86400 * 24
  lat = 40.0 * pi / 180  # Latitude (degrees -> radians) for solar radiation

  ## Meteorological forcing 
  # Ta = 20.0
  Ta = gen_Ta(hour)     # Air temperature (C)
  es, d_es = satvap(Ta) # [degC] -> Pa
  ea = es * RH / 100
  update_met!(met, Ta, ea, u, rain, snow, Pa)

  ## Flux
  Ts = Tmean + K0 # [degK]
  es, d_es = satvap(Ts - K0)
  flux = Flux{Float64}(; θ_surf=Ts, e_surf=es) # 生成一次即可

  ## Radiation and Canopy
  rad = Radiation()
  update_rad!(rad, day, hour, lat, Ta)
  # rad = Radiation(day, hour, lat; Rln_in=_cal_Rli(Ta))

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
