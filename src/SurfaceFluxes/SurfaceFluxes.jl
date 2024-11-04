module SurfaceFluxes

using HydroTools
using HydroTools: Γd, λ_fus, λ_vap, λ_sub
using HydroTools: root_hybrid, root_brent
# import HydroTools: soil_temperature_delta

using Parameters: @with_kw
using DocStringExtensions
using UnPack

export soil_depth_init
export soil_temperature, soil_temperature_delta, soil_thermal_properties
export ρ_wat, ρ_ice

export Met, Radiation, Canopy, Soil, Flux
export cal_coszen
export satvap
export Soil, init_soil!
export MOST, surface_fluxes!

# λ_fus = 0.3337e6                     # Heat of fusion for water at 0 C (J/kg)
ρ_wat = 1000.0                       # Density of water (kg/m3)
ρ_ice = 917.0                        # Density of ice (kg/m3)
tfrz = K0                           # freezing Temperature (K)

M_h2o = 18.02 * 1e-3 # [kg mol-1]
M_dry = 28.97 * 1e-3

module physcon
# vkc    = 0.4             # von Karman constant
# grav   = 9.80665         # Gravitational acceleration (m/s2)
tfrz = 273.15          # Freezing point of water (K)
σ = 5.67e-08        # Stefan-Boltzmann constant (W/m2/K4)

ρ_wat = 1000.0          # Density of water (kg/m3)
ρ_ice = 917.0           # Density of ice (kg/m3)
cp_wat = 4188.0          # Specific heat of water (J/kg/K)
cp_ice = 2117.27         # Specific heat ice (J/kg/K)
cv_wat = cp_wat * ρ_wat  # Heat capacity of water (J/m3/K)
cv_ice = cp_ice * ρ_ice  # Heat capacity of ice (J/m3/K)
tk_wat = 0.57            # Thermal conductivity of water (W/m/K)
tk_ice = 2.29            # Thermal conductivity of ice (W/m/K)
end

include("Soil.jl")
include("soil_depth_init.jl")
include("soil_thermal_properties.jl")
include("soil_thermal_properties_flux.jl")
include("soil_temperature.jl")
include("soil_temperature_delta.jl")

# using ProtoStructs
include("Met.jl")
include("Radiation.jl")
include("Canopy.jl")
include("Flux.jl")
# 三个高度：ref, surf, ground：参考高度，0表面，
include("tools_met.jl")
include("MOST.jl")
include("surface_fluxes.jl")

end
