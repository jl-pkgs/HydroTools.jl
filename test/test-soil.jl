@testset "soil" begin
  dz = [0.1, 0.2, 0.3]
  nsoil, z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  @test z ≈ [-0.05, -0.2, -0.45]
  @test z₊ₕ ≈ -cumsum(dz)
  # soil_temperature(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
end
