@testset "soil" begin
  dz = [0.1, 0.2, 0.3]
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  @test z ≈ [-0.05, -0.2, -0.45]
  @test z₊ₕ ≈ -cumsum(dz)
  # soil_temperature(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
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

  m_sat = Θ_S[k] * ρ_wat * dz # kg/m2
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
    Tsoil_next2, G2 = soil_temperature(dz, dt, κ, cv, Tsoil_cur, Tsurf_next; method, solution = "crank-nicolson")
    @test G ≈ 464.6218567684632
    @test G2 ≈ 392.1533026315689
  end
end
