
"""
    Cambell(ψ; param)

Campbell (1974) relationships

# Arguments
+ `ψ`: Matric potential, cm
+ `param`
  - `ψ_sat`: Matric potential at saturation, [cm]
  - `θ_sat`: Volumetric water content at saturation
  - `b`    : Exponent
  - `K_sat`: Hydraulic conductivity at saturation [cm h-1]

# TODO: 核对变量的单位

# Examples
```julia
param = (θ_sat = 0.25, ψ_sat = -25.0, b = 0.2, K_sat = 3.4e-03)
θ, K, ∂θ∂ψ = Cambell(-100; param)
```
"""
@fastmath function Cambell(ψ::Real; param)
  @unpack ψ_sat, θ_sat, b, K_sat = param

  # Volumetric soil moisture (θ) for specified matric potential 
  θ = ψ < ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : ψ_sat

  # Hydraulic conductivity (K) for specified matric potential
  K = ψ < ψ_sat ? K_sat * (θ / θ_sat)^(2b + 3) : K_sat

  # Cap = dθ/dψ
  ∂θ∂ψ = ψ < ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : ψ_sat

  θ, K, ∂θ∂ψ
end

