using Test
using HydroTools


include("test-thermal.jl")
include("test-soil_moisture.jl")
include("test-soil_temperature.jl")
include("test-radiation.jl")
include("test-PMLV2.jl")
include("test-sceua.jl")

@testset "HW_index" begin
  anorm = [0.1, 0.2, 0.3, 0.2, 0.1, 0, -0.1, 0.1, 0.2, 0.3]
  res = HW_index(anorm)
  @test res.duration == 9
  @test res.frequency == 2
  @test res.intensity == 0.3
  @test res.volume == 1.5
  @test res.PR ≈ 89.99999999999993
  @test res.FAR ≈ 0.9888888888888889

  @test HW_index([0.1, 0.2, 0.3, 0.2, 0.1, 0, -0.1, 0.1, 0.2, 0.3]) ==
        (duration=9, frequency=2, intensity=0.3, volume=1.5, PR=89.99999999999993, FAR=0.9888888888888889)
  # @test HW_index([-1, -1]) ==
  #       (duration=0, frequency=0, intensity=NaN, volume=NaN, PR=NaN, FAR=NaN)
end

@testset "heat_index" begin
  @test heat_index(30., 50.) == 31.049081444444305
  @test heat_index(28.3, 87.0) == 34.22188833448005
  @test heat_index(28.3, 12.0) == 26.761397287925213
  @test heat_index(30., 40.) < heat_index(30., 70.)
end

@testset "ET0 models" begin
  @test ET0_eq(200.0, 20.0, 2.0) == (2.4536, 0.14474018811241365, 0.0013275295705686334, 6.978705385597235)
  @test ET0_Penman48(200.0, 20.0, 2.0, 2.0) ≈ 8.253727734102089
  @test ET0_FAO98(200.0, 20.0, 2.0, 2.0) ≈ 7.159788896344576
end


@testset "detect_events" begin
  y = [0, 1, 2, 3, 0, 2, 1, 0]
  lgl = y .> 0
  @test detect_events(lgl) == [(i_beg=2, i_end=4, len=3), (i_beg=6, i_end=7, len=2)]

  @test detect_events(y, lgl) == [
    (i_beg=2, i_peak=4, i_end=4, len=3, len_left=2, len_right=0, peak=3)
    (i_beg=6, i_peak=6, i_end=7, len=2, len_left=0, len_right=1, peak=2)
  ]
end

@testset "GOF" begin
  @test GOF(1:10, 2:11) ==
        (NSE=0.8787878787878788, R2=1.0, KGE=0.8181818181818181, R=1.0, RMSE=1.0, MAE=1.0, bias=1.0, bias_perc=18.181818181818183, n_valid=10)
end
