@testset "ET PMLV2" begin
  Prcp = 2.0 # mm
  Tavg = 20.0
  Rs = 200.0
  Rn = 50.0
  VPD = 2.0
  U2 = 2.0
  LAI = 2.0

  par = param0
  par = add(par, list(hc=2.0))

  @test_nowarn PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI; par)
end
