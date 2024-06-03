using HydroTools
using HydroTools.SurfaceFluxes
using Test

begin
  z = 30.0 # Reference height (m)
  hc = 20.0 # canopy height
  d = 0.67 * hc # Zero-plane displacement height (m)
  z0m = 0.13 * hc # Roughness length for momentum (m)
  z0c = 0.10 * z0m # Roughness length for scalars (m)
  param = (; z, d, z0m, z0c)
end

day = 1
hour = 10
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
flux = Flux(; θ_surf=Ts, e_surf=es)

@testset "MOST" begin
  @test MOST(10.0, met, flux; param...) ≈ -10.745758995719331  
end
