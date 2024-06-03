module SurfaceFluxes

using HydroTools
using HydroTools: Γd, λ_fus, λ_vap, λ_sub
using HydroTools: root_hybrid, root_brent
import HydroTools: soil_temperature_delta

using Parameters: @with_kw
using DocStringExtensions
using UnPack


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

include("DataType/DataType.jl")
include("tools_met.jl")
include("MOST.jl")
include("surface_fluxes.jl")

export satvap
export Soil, init_soil!
export MOST, surface_fluxes

end
