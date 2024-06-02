using HydroTools
using Test

@testset "thermal" begin
  @test theta_wet(850.0, 20.0, 18.0).Θ_se ≈ 76.23036568899175
  @test theta_wet_bolton(850.0, 20.0, 18.0).Θ_se ≈ 80.67244668155882
  # @test_nowarn cal_Rn(Rs, Rln, Tavg, albedo, emiss)
end
