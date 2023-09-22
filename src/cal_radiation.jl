"""
    get_Rn(Rs, Rln, Tavg, albedo, emiss)

Calculate net radiation (Rn) using input parameters.

# Arguments

- `Rs`   : Shortwave inward solar radiation, W m-2.
- `Rln`  : Longwave inward solar radiation, W m-2.
- `Tavg` : Average temperature in Celsius.
- `α`    : Shortwave albedo.
- `ϵ`    : Longwave emissivity.

# Returns
- `Rn   `: Net radiation, W m-2.
"""
function cal_Rn(Rs::Quantity{T}, Rln::Quantity{T}, Tavg::Quantity{T},
  α::Float64, ϵ::Float64) where {T<:Real}
  
  K = uconvert(u"K", Tavg)
  RLout = Stefan * K^4 |> unit(Rs)

  Rnl = ϵ * (Rln - RLout)
  Rns = (1.0 - α) * Rs
  Rn = Rns + Rnl

  max(Rn, 0.0 * unit(Rn))
end

# 使用Unit速度会有较大牺牲
function cal_Rn(Rs::T, Rln::T, Tavg::T, α::Float64, ϵ::Float64) where {T<:Real}
  RLout = Stefan.val * (Tavg + 273.15)^4 / 0.0864 # convert from MJ m-2 d-1 to W/m2

  Rnl = ϵ * (Rln - RLout)
  Rns = (1.0 - α) * Rs
  Rn = Rns + Rnl
  max(Rn, 0.0)
end


export cal_Rn
