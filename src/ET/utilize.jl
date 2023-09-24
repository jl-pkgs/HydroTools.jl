function PMLV2(Prcp::Vector{T}, Tavg::Vector{T},
  Rs::Vector{T}, Rn::Vector{T},
  VPD::Vector{T}, U2::Vector{T}, LAI::Vector{T},
  Pa::Vector{T}, Ca=380.0;
  par=param0, frame=3,
  res::Union{Nothing,output_PML}=nothing) where {T<:Real}

  n = length(Prcp)
  fields = fieldnames(interm_PML)[2:end-2]
  res === nothing && (res = interm_PML{T}(; n))

  for t = 1:n
    r = PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t], Ca; par)
    res[t, fields] = r
  end
  
  fval_soil = movmean2(Prcp, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(fval_soil, 0, 1)

  res.Es .= fval_soil .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end

# if option == 1
#   res.GPP[t], res.Ec[t], res.Ecr[t], res.Eca[t],
#   res.Ei[t], res.Pi[t], res.Es_eq[t], res.Eeq[t], res.ET_water[t],
#   res.Ga[t], res.Gc_w[t] =
#     PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t], Ca; par)
# elseif option == 2
# end
