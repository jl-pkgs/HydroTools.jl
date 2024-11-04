import HydroTools.SurfaceFluxes: θ_S

@testset "soil" begin
  dz = [0.1, 0.2, 0.3]
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  @test z ≈ [-0.05, -0.2, -0.45]
  @test z₊ₕ ≈ -cumsum(dz)
  # soil_temperature(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
end


@testset "soil_temperature_delta" begin
  dz = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058]
  dt = 1800 # seconds, 0.5h
  n = 10
  κ = fill(1.4651477226706402, n)
  cv = fill(2.6628438e6, n)
  Tsoil_cur = fill(K0 + 25, n)
  df0 = -148.3184062187158
  f0 = -798.1091814317192

  Tsoil_next, G_soil, G_snow = soil_temperature_delta(dz, dt, κ, cv, Tsoil_cur, df0, f0)
  @test Tsoil_next[1] ≈ 294.30394328396324
  @test Tsoil_next[2] ≈ 296.2755002132709
  @test Tsoil_next[3] ≈ 297.56178599804434
end


# for itime = 1:ntime
@testset "soil_temperature" begin
  ## Constants ===================================================================
  tmean = K0 + 15
  trange = 10

  dt = 1800 # seconds, 0.5h
  ntime = round(86400 / dt)
  nday = 200

  soil_texture = 1
  k = soil_texture
  ## =============================================================================

  n = 120
  dz = fill(0.025, n)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  m_sat = θ_S[k] * ρ_wat * dz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat

  begin
    Tsoil = fill(2.0 + K0, n)

    itime = 1
    hour = itime * (dt / 86400 * 24)
    Tsurf_next = tmean + 0.5 * trange * sin(2 * pi / 24 * (hour - 8.0)) # 采用正弦函数来反映温度变化
    Tsoil_cur = deepcopy(Tsoil)

    method = "apparent-heat-capacity"
    κ, cv = soil_thermal_properties(dz, Tsoil, m_liq, m_ice; soil_texture=1, method)
    Tsoil_next, G = soil_temperature(dz, dt, κ, cv, Tsoil_cur, Tsurf_next; method)
    @test G ≈ 464.6218567684632

    Tsoil_next2, G2 = soil_temperature(dz, dt, κ, cv, Tsoil_cur, Tsurf_next; method, solution = "crank-nicolson")
    # 已核对，与MATLAB的版本结果一致
    @test G2 ≈ 392.1533026315689
  end
  
  ## 这里需要做一个测试，土壤温度多久可以达到稳态
end
