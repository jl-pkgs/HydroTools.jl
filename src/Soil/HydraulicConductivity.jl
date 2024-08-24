using UnPack

include("ψ_Cambell.jl")
include("ψ_van_Genuchten.jl")

export Cambell, van_Genuchten
export matric_potential, HydraulicConductivity

# TODO: add a parameter struct for `matric_potential`

# !deprecated
function HydraulicConductivity(ψ::AbstractVector; param, fun=van_Genuchten)
  n = length(ψ)
  Θ = zeros(n)
  K = zeros(n)
  ∂Θ∂ψ = zeros(n)
  for i in 1:n
    Θ[i], K[i], ∂Θ∂ψ[i] = fun(ψ[i]; param)
  end
  Θ, K, ∂Θ∂ψ
end


# Calculate ψ for a given Θ
function matric_potential(Θ, param; method="van_Genuchten")
  if method == "van_Genuchten"
    @unpack Θ_res, Θ_sat, α, n, m = param
    Se = @. (Θ - Θ_res) / (Θ_sat - Θ_res)
    ψ = @. -((Se^(-1 / m) - 1)^(1 / n)) / α
  elseif method == "Campbell"
    @unpack ψ_sat, Θ_sat, b = param
    ψ = @. ψ_sat * (Θ / Θ_sat)^-b
  end
  ψ
end
