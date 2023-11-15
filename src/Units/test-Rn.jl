## 使用Unit速度会有较大牺牲
# function cal_Rn(Rs::Quantity{T}, Rln::Quantity{T}, Tavg::Quantity{T},
#   α::Float64, ϵ::Float64) where {T<:Real}
#   K = uconvert(u"K", Tavg)
#   RLout = Stefan * K^4 |> unit(Rs)
#   Rnl = ϵ * (Rln - RLout)
#   Rns = (1.0 - α) * Rs
#   Rn = Rns + Rnl
#   max(Rn, 0.0 * unit(Rn))
# end
