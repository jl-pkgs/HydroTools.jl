"""
    van_Genuchten(ψ, param)

van Genuchten (1980) relationships

# Arguments
+ `ψ`: Matric potential
+ `param`
  - `θ_res`       : Residual water content
  - `θ_sat`       : Volumetric water content at saturation
  - `α`           : Inverse of the air entry potential (cm-1)
  - `n`           : Pore-size distribution index
  - `m`           : Exponent
  - `K_sat`       : Hydraulic conductivity at saturation (cm/s)
  - `soil_texture`: Soil texture flag

# Examples
```julia
# Haverkamp et al. (1977): sand
param = (soil_texture = 1, 
  θ_res = 0.075, θ_sat = 0.287, 
  α = 0.027, n = 3.96, m = 1, K_sat = 34 / 3600)

# Haverkamp et al. (1977): Yolo light clay
param = (soil_texture=2, 
  θ_res = 0.124, θ_sat = 0.495,
  α = 0.026, n = 1.43, m = 1 - 1 / 1.43,
  K_sat = 0.0443 / 3600)
```
"""
function van_Genuchten(ψ::Real; param)
  @unpack θ_res, θ_sat, α, n, m, K_sat, soil_texture = param

  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1

  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  θ = θ_res + (θ_sat - θ_res) * Se

  # Hydraulic conductivity (K) for specified matric potential (ψ)
  if Se <= 1
    K = K_sat * sqrt(Se) * (1 - (1 - Se^(1 / m))^m)^2
    # Special case for Haverkamp et al. (1977) sand (soil_texture = 1) and Yolo light clay (soil_texture = 2)
    if soil_texture == 1
      K = K_sat * 1.175e6 / (1.175e6 + abs(ψ)^4.74)
    elseif soil_texture == 2
      K = K_sat * 124.6 / (124.6 + abs(ψ)^1.77)
    end
  else
    K = K_sat
  end

  # Specific moisture capacity (∂θ∂ψ) for specified matric potential (ψ)
  if ψ <= 0.0
    num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n - 1)
    den = (1 + (α * abs(ψ))^n)^(m + 1)
    ∂θ∂ψ = num / den
  else
    ∂θ∂ψ = 0.0
  end
  θ, K, ∂θ∂ψ
end
