using HydroTools
using Statistics
using Random
using Test


@testset "index_prcp" begin
  Random.seed!(1)
  x = rand(1000) * 10
  # q95 = quantile(x, 0.95)
  r = index_prcp(x) # cdd, r20mm, rx5day, r95ptot, prcptot
  
  @test round.(r, digits=2) == [2.0, 0, 10.0, 489.80, 4939.32] 
end


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
  @test heat_index(30.0, 50.0) == 31.049081444444305
  @test heat_index(28.3, 87.0) == 34.22188833448005
  @test heat_index(28.3, 12.0) == 26.761397287925213
  @test heat_index(30.0, 40.0) < heat_index(30.0, 70.0)
end
