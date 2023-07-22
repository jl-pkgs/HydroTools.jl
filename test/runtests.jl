using Test
using HydroTools


@testset "HW_index" begin
  anorm = [0.1, 0.2, 0.3, 0.2, 0.1, 0, -0.1, 0.1, 0.2, 0.3]
  res = HW_index(anorm)
  @test res.duration == 9
  @test res.frequency == 2
  @test res.intensity == 0.3
  @test res.volume == 1.5
  @test res.PR ≈ 89.99999999999993
  @test res.FAR ≈ 0.9888888888888889
end

