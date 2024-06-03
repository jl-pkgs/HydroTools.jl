"""
# Fields
$(TYPEDFIELDS)
"""
@with_kw mutable struct Flux{T<:AbstractFloat}
  "temperature of surface (K)"
  θ_surf::T = 0.0
  "vapor pressure of surface (Pa)"
  e_surf::T = 0.0

  u₊::T = 0.0
  θ₊::T = 0.0
  q₊::T = 0.0

  "Obukhov length"
  ζ::T = 0.0
  "Aerodynamic conductances, [mol/m2/s]"
  g_ac::T = 0.0

  "Latent heat flux (W m-2)"
  LE::T = 0.0
  "Sensible heat flux (W m-2)"
  H::T = 0.0
  "Ground Surface Soil heat flux (W m-2)"
  G_soil::T = 0.0
  "Ground Surface snow heat flux (W m-2)"
  G_snow::T = 0.0
end
