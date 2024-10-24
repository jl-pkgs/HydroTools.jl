using HydroTools, Test


@testset "face and center" begin
  dz = [2, 4, 8, 10]
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  @test F2C(z₊ₕ) == z
  @test C2F(z) == z₊ₕ
end

@testset "Cambell" begin
  param = (θ_sat=0.25, ψ_sat=-25.0, b=0.2, K_sat=3.4e-03)
  θ, K, ∂θ∂ψ = Cambell(-100; param)
  @test (θ, K, ∂θ∂ψ) == (0.000244140625, 1.979060471057893e-13, 1.220703125e-5)
end

@testset "soil_moisture!" begin
  n = 150
  dz = ones(n)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  soil_texture = 1
  param = (soil_texture=1,
    θ_res=0.075, θ_sat=0.287,
    α=0.027, n=3.96, m=1, K_sat=34 / 3600)

  θ = fill(0.1, n)
  ψ = matric_potential(θ, param; method="van_Genuchten")

  θ0 = 0.267
  ψ0 = matric_potential(θ0, param; method="van_Genuchten")

  dt = 5
  ntim = 0.8 * 3600 / dt

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  @time for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\n", hour)
    
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture!(θ, ψ, ψ0, dz, dt, param)

    # % Sum fluxes for relative mass balance error
    sum_in += abs(Q0) * dt
    sum_out += abs(QN) * dt
    sum_store += dθ
  end

  @test sum_in ≈ 11.810243822643141
  @test sum_out ≈ 0.10508872215771699
  @test sum_store ≈ 11.704825251924781
end
