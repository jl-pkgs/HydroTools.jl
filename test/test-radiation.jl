using Test

@testset "radiation" begin
  cal_Rsi()
  cal_Rsi_toa()

  cal_Rln(30, 20, 1, 0.5)
  cal_Rln_out(10)
  cal_Rli(10)
  cal_Rln_yang2019(20, 100, 200)


  Rs = 200
  Rln = 400
  Tavg = 20
  albedo = 0.2
  emiss = 0.96
  @test_nowarn cal_Rn(Rs, Rln, Tavg, albedo, emiss)
end
