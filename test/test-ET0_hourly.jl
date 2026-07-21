@testset "Hourly ET0" begin
  hour = Float32.(0:23)
  Ta = 20f0 .+ 6f0 .* sinpi.((hour .- 9f0) ./ 12f0)
  Tdew = 12f0 .+ 2f0 .* sinpi.((hour .- 9f0) ./ 12f0)
  Rn_G = max.(0f0, sinpi.((hour .- 6f0) ./ 12f0))
  U2 = 2f0 .+ 0.5f0 .* sinpi.((hour .- 6f0) ./ 12f0) .^ 2
  Pa = fill(101.325f0, 24)
  VPD = max.(cal_es.(Ta) .- cal_es.(Tdew), 0f0)

  fao = HydroTools._FAO98.(Rn_G, Ta, VPD, U2, Pa)
  pt = HydroTools._PT1972.(Rn_G, Ta, Pa)
  romanenko = HydroTools._Romanenko.(Ta, Tdew)

  @test eltype(fao) == eltype(pt) == eltype(romanenko) == Float32
  @test all(isfinite, fao) && all(>=(0), fao)
  @test all(isfinite, pt) && all(>=(0), pt)
  @test all(isfinite, romanenko) && all(>=(0), romanenko)
  @test sum(fao) ≈ 3.45895f0 rtol=1f-4
  @test sum(pt) ≈ 2.81174f0 rtol=1f-4
  @test sum(romanenko) ≈ 6.00723f0 rtol=1f-4
  @test all(iszero, pt[Rn_G .== 0])
  @test HydroTools._Romanenko(20f0, 20f0) == 0f0
end
