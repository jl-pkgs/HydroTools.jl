using UnPack

include("ψ_Cambell.jl")
include("ψ_van_Genuchten.jl")

export Cambell, van_Genuchten
export matric_potential, HydraulicConductivity

# TODO: add a parameter struct for `matric_potential`

# !deprecated
function HydraulicConductivity(ψ::AbstractVector; param, fun=van_Genuchten)
  n = length(ψ)
  θ = zeros(n)
  K = zeros(n)
  ∂θ∂ψ = zeros(n)
  for i in 1:n
    θ[i], K[i], ∂θ∂ψ[i] = fun(ψ[i]; param)
  end
  θ, K, ∂θ∂ψ
end


# Calculate ψ for a given θ
function matric_potential(θ, param; method="van_Genuchten")
  if method == "van_Genuchten"
    @unpack θ_res, θ_sat, α, n, m = param
    Se = @. (θ - θ_res) / (θ_sat - θ_res)
    ψ = @. -((Se^(-1 / m) - 1)^(1 / n)) / α
  elseif method == "Campbell"
    @unpack ψ_sat, θ_sat, b = param
    ψ = @. ψ_sat * (θ / θ_sat)^-b
  end
  ψ
end
