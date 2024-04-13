# using MakieLayers
# using GLMakie

function _cal_PET!(
  lat::AbstractVector, doys::AbstractVector,
  Tmin::V3, Tmax::V3, Tavg::V3, RH::V3, SSD::V3,
  U2::V3, alt::V2) where {
  V3<:AbstractArray{<:Real,3},V2<:AbstractArray{<:Real,2}}

  PET = zeros(Float32, size(Tmin)...) .* NaN32
  Rn = zeros(Float32, size(Tmin)...) .* NaN32

  @inbounds @showprogress for i = eachindex(lon)
    for j = eachindex(lat)
      Z = alt[i, j]
      pa = cal_Pa(Z) |> Float32
      _lat = lat[j] |> Float32

      @views begin
        u2 = U2[i, j, :]
        rh = RH[i, j, :]
        tavg = Tavg[i, j, :]
        tmin = Tmin[i, j, :]
        tmax = Tmax[i, j, :]
        ssd = SSD[i, j, :]
      end
      all(isnan(tavg)) && continue

      ea = cal_ea.(tavg, rh)
      vpd = cal_es.(tavg) .- ea

      # cal_Rn(lat, J, Tmin::T, Tmax::T, ea::T, ssd::T; albedo, Z)
      _Rn = cal_Rn.(_lat, doys, tmin, tmax, ea, ssd; Z)
      _pet = ET0_FAO98.(_Rn, tavg, vpd, u2, pa) # mm/d

      Rn[i, j, :] .= _Rn
      PET[i, j, :] .= _pet
    end
  end
  PET, Rn
end
