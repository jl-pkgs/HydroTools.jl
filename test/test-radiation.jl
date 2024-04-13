using Test

@testset "radiation" begin
  cal_Rsi(20.0, 20, 8.0)
  cal_Rsi_toa()

  cal_Rnl(30., 20., 1.0, 0.5)
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


@testset "Norman_Longwave" begin
  n = 50
  ϵ = [1.0; fill(0.98, n - 1)]
  T_leaf = [20.0; fill(25.0, n - 1)]
  τd = ones(n) .* 0.915     # transmittance of diffuse radiation through each leaf layer
  L_up, L_dn, Rln, Rln_soil, Rln_veg = Norman_Longwave(T_leaf, ϵ, τd)
  @test Rln_soil ≈ 28.382474037679856
  @test Rln_veg ≈ -75.5489070386997
end

@testset "Norman_Shortwave" begin
  PAR_sun, PAR_sha, frac_sha, frac_sun = Norman_Shortwave([0.1, 0.2, 0.3])
  @test PAR_sha ≈ [168.89332828782798, 159.20858966542008, 143.23329123213492]
end

@testset "cal_Rn" begin
  r = cal_Rn(20.0, 20, 20.0, 25.0, 2.0, 10.0) # default unit is W m-2
  @test W2MJ(r.Rn) ≈ 9.91692555
end
