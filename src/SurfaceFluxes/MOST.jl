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
  (; θ, θv, e, u, Pa, M_air) = met
  (; θ_surf, e_surf) = flux
  ζ0 = abs(ζ0) <= 0.1 ? 0.1 : ζ0 # Prevent near-zero values
  z_minus_d = z - d

  # Evaluate psi for momentum at the reference height (zref-disp) and surface (z0m)
  Ψm = -Ψ_m_monin_obukhov(z0m / ζ0, z_minus_d / ζ0) # Eq 6.39
  Ψc = -Ψ_c_monin_obukhov(z0c / ζ0, z_minus_d / ζ0) # Eq 6.39

  zlog_m = log(z_minus_d / z0m)
  zlog_c = log(z_minus_d / z0c)

  # Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)
  u₊ = u * k / (zlog_m + Ψm)                          # Eq. 6.40
  θ₊ = (θ - θ_surf) * k / (zlog_c + Ψc)               # Eq. 6.43
  q₊ = (e - e_surf) / Pa * k / (zlog_c + Ψc)          # Eq. 6.44
  θv₊ = θ₊ + 0.61 * θ * q₊ * (M_h2o / M_air)          # Eq. 6.33

  ζ = u₊^2 * θv / (k * g * θv₊)                        # Eq. 6.30
  @pack! flux = u₊, θ₊, q₊, ζ
  return ζ0 - ζ # dζ
end

# ϕ: phi
# Ψ: psi
function phi_m_monin_obukhov(x)
  # Evaluate the Monin-Obukhov phi function for momentum at x
  # Eq. 6.37
  if x < 0
    phi_m = (1 - 16 * x)^(-0.25)
  else
    phi_m = 1 + 5 * x
  end
  return phi_m
end

function phi_c_monin_obukhov(x::Real)
  # Evaluate the Monin-Obukhov phi function for scalars at x
  # Eq. 6.38
  if x < 0
    phi_c = (1 - 16 * x)^(-0.5)
  else
    phi_c = 1 + 5 * x
  end
  return phi_c
end


function Ψ_c_monin_obukhov(x::Real)
  # Evaluate the Monin-Obukhov Ψ function for scalars at x, Eq. 6.47
  if x < 0
    y = (1 - 16 * x)^0.25
    Ψ_c = 2 * log((1 + y^2) / 2)
  else
    Ψ_c = -5 * x
  end
  return Ψ_c
end

function Ψ_m_monin_obukhov(x)
  # Evaluate the Monin-Obukhov Ψ function for momentum at x, Eq. 6.46
  if x < 0
    y = (1 - 16 * x)^0.25
    Ψ_m = 2 * log((1 + y) / 2) + log((1 + y^2) / 2) - 2 * atan(y) + pi / 2
  else
    Ψ_m = -5 * x
  end
  return Ψ_m
end


function Ψ_c_monin_obukhov(x1::T, x2::T) where {T<:Real}
  Ψ_c_monin_obukhov(x2) - Ψ_c_monin_obukhov(x1)
end

function Ψ_m_monin_obukhov(x1::T, x2::T) where {T<:Real}
  Ψ_m_monin_obukhov(x2) - Ψ_m_monin_obukhov(x1)
end
