"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw struct Soil{T<:AbstractFloat}
  "layers of soil"
  n::Int = 5
  "time step (s) of solver"
  dt::T = 1800.0
  "depth of soil layer (m)"
  dz::Vector{T} = zeros(T, n)
  "Soil depth (m) at center of layer i (negative distance from surface)"
  z::Vector{T} = zeros(T, n)
  "Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface), Fig. 5.3"
  z₊ₕ::Vector{T} = zeros(T, n)
  "Thickness between between z(i) and z(i+1)"
  dz₊ₕ::Vector{T} = zeros(T, n)

  "temperature of soil layer"
  Tsoil::Vector{T} = zeros(T, n)

  "heat capacity, [J m-3 K-1]"
  cv::Vector{T} = zeros(T, n)
  "thermal conductivity, [W m-1 K-1]"
  κ::Vector{T} = zeros(T, n)

  "soil moisture of ice, [kg m-2]"
  SM_ice::Vector{T} = zeros(T, n)
  "soil moisture of liquid, [kg m-2]"
  SM_liq::Vector{T} = zeros(T, n)
end

function Soil(dz; dt=1800, kw...)
  n = length(dz)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  Soil{Float64}(; n, dt, dz, z, z₊ₕ, dz₊ₕ, kw...)
end


function init_soil!(soil::Soil; Ts=25.0, soil_texture::Int=5, satfrac=0.85)
  (; n, κ, cv, Tsoil, SM_liq, SM_ice, dz) = soil
  ρ_wat = 1000.0          # Density of water (kg/m3)
  isa(Ts, Number) && (Ts = fill(Ts, n))
  
  # Initial soil temperature (K) and unfrozen and frozen water (kg H2O/m2)
  for i in 1:n
    Tsoil[i] = Ts[i] + K0 # Temperature
    # Soil water at saturation (kg H2O/m2)
    SM_sat = θ_S[soil_texture] * ρ_wat * dz[i]
    # Actual water content is some fraction of saturation. These are only used for soil
    # thermal properties and phase change. Note the inconsistency with the use of soil
    # water in the bucket model to calculate the soil wetness factor.
    if Tsoil[i] > K0
      SM_ice[i] = 0
      SM_liq[i] = satfrac * SM_sat
    else
      SM_ice[i] = satfrac * SM_sat
      SM_liq[i] = 0
    end
  end
  soil_thermal_properties_flux(κ, cv, Tsoil, SM_liq, SM_ice, dz; soil_texture)
  soil
end


function soil_temperature_delta(soil::Soil, df0::Real, f0::Real, snow_water::Real=0.0)
  (; dz, dt, κ, cv, Tsoil) = soil
  soil_temperature_delta(dz, dt, κ, cv, Tsoil, df0, f0, snow_water) # Tsoil_next, G_soil, G_snow
end
