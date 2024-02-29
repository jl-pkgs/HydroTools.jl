using UnPack

"""
    Cambell(ϕ; param)

Campbell (1974) relationships

# Arguments
+ `ϕ`: Matric potential
+ `param`
  - `ϕ_sat`: Matric potential at saturation
  - `Θ_sat`: Volumetric water content at saturation
  - `b`    : Exponent
  - `K_sat`: Hydraulic conductivity at saturation

# Examples
```julia
param = (Θ_sat = 0.25, ϕ_sat = -25.0, b = 0.2, K_sat = 3.4e-03)
```
"""
@fastmath function Cambell(ϕ::Real; param)
  @unpack ϕ_sat, Θ_sat, b, K_sat = param

  # Volumetric soil moisture (Θ) for specified matric potential 
  Θ = ϕ < ϕ_sat ? Θ_sat * (ϕ / ϕ_sat)^(-1 / b) : ϕ_sat

  # Hydraulic conductivity (K) for specified matric potential
  K = ϕ < ϕ_sat ? K_sat * (Θ / Θ_sat)^(2b + 3) : K_sat

  # Cap = dΘ/dϕ = 
  ∂Θ∂ϕ = ϕ < ϕ_sat ? -Θ_sat / (b * ϕ_sat) * (ϕ / ϕ_sat)^(-1 / b - 1) : ϕ_sat

  Θ, K, ∂Θ∂ϕ
end


"""
    van_Genuchten(ϕ, param)

van Genuchten (1980) relationships

# Arguments
+ `ϕ`: Matric potential
+ `param`
  - `Θ_res`       : Residual water content
  - `Θ_sat`       : Volumetric water content at saturation
  - `α`           : Inverse of the air entry potential (cm-1)
  - `n`           : Pore-size distribution index
  - `m`           : Exponent
  - `K_sat`       : Hydraulic conductivity at saturation (cm/s)
  - `soil_texture`: Soil texture flag

# Examples
```julia
# Haverkamp et al. (1977): sand
param = (soil_texture = 1, 
  Θ_res = 0.075, Θ_sat = 0.287, 
  α = 0.027, n = 3.96, m = 1, K_sat = 34 / 3600)

# Haverkamp et al. (1977): Yolo light clay
param = (soil_texture=2, 
  Θ_res = 0.124, Θ_sat = 0.495,
  α = 0.026, n = 1.43, m = 1 - 1 / 1.43,
  K_sat = 0.0443 / 3600)
```
"""
function van_Genuchten(ϕ::Real; param)
  @unpack Θ_res, Θ_sat, α, n, m, K_sat, soil_texture = param

  # Effective saturation (Se) for specified matric potential (ϕ)
  Se = ϕ <= 0 ? (1 + (α * abs(ϕ))^n)^-m : 1

  # Volumetric soil moisture (Θ) for specified matric potential (ϕ)
  Θ = Θ_res + (Θ_sat - Θ_res) * Se

  # Hydraulic conductivity (K) for specified matric potential (ϕ)
  if Se <= 1
    K = K_sat * sqrt(Se) * (1 - (1 - Se^(1 / m))^m)^2
    # Special case for Haverkamp et al. (1977) sand (soil_texture = 1) and Yolo light clay (soil_texture = 2)
    if soil_texture == 1
      K = K_sat * 1.175e6 / (1.175e6 + abs(ϕ)^4.74)
    elseif soil_texture == 2
      K = K_sat * 124.6 / (124.6 + abs(ϕ)^1.77)
    end
  else
    K = K_sat
  end

  # Specific moisture capacity (∂Θ∂ϕ) for specified matric potential (ϕ)
  if ϕ <= 0.0
    num = α * m * n * (Θ_sat - Θ_res) * (α * abs(ϕ))^(n - 1)
    den = (1 + (α * abs(ϕ))^n)^(m + 1)
    ∂Θ∂ϕ = num / den
  else
    ∂Θ∂ϕ = 0.0
  end
  Θ, K, ∂Θ∂ϕ
end


function HydraulicConductivity(ϕ::AbstractVector; param, fun=van_Genuchten)
  n = length(ϕ)
  Θ = zeros(n)
  K = zeros(n)
  ∂Θ∂ϕ = zeros(n)
  for i in 1:n
    Θ[i], K[i], ∂Θ∂ϕ[i] = fun(ϕ[i]; param)
  end
  Θ, K, ∂Θ∂ϕ
end


# Calculate ϕ for a given Θ
function matric_potential(Θ, param; method="van_Genuchten")
  if method == "van_Genuchten"
    @unpack Θ_res, Θ_sat, α, n, m = param
    Se = @. (Θ - Θ_res) / (Θ_sat - Θ_res)
    ϕ = @. -((Se^(-1 / m) - 1)^(1 / n)) / α
  elseif method == "Campbell"
    @unpack ϕ_sat, Θ_sat, b = param
    ϕ = @. ϕ_sat * (Θ / Θ_sat)^-b
  end
  ϕ
end


export matric_potential
