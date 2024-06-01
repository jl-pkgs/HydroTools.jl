"""
    soil_temperature_delta(dz, dt, κ, cv, Tsoil_cur, Tsurf_next)

> 这里采用的是方案1的边界条件, sum(F) = G, G = - k1 * (T1 - T0) / dz1

# Arguments: 

- `κ`: thermal conductivity (W/m/K)
- `cv`: volumetric heat capacity (J/m3/K)
- `f0`: 
- `Tsurf_next`: Tsurf_next_next, T0_{n+1}
- `solution`:
  + `implicit`:
  + `crank-nicolson`:

# Examples

```julia
dz = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058]
dt = 1800 # seconds, 0.5h
n = 10
κ = fill(1.4651477226706402, n)
cv = fill(2.6628438e6, n)
Tsoil_cur = fill(K0 + 25, n)
df0 = -148.3184062187158
f0 = -798.1091814317192

Tsoil_next, G = soil_temperature_delta(dz, dt, κ, cv, Tsoil_cur, df0, f0, snow_water)
```
"""
function soil_temperature_delta(dz::AbstractVector, dt::Real, 
  κ::AbstractVector, cv::AbstractVector, Tsoil_cur::AbstractVector,
  df0::Real, f0::Real, snow_water::Real=0.0)
  # solution = "implicit", method = "apparent-heat-capacity"
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  n = length(dz)

  # Thermal conductivity at interface (W/m/K)
  κ₊ₕ = zeros(1, n - 1)
  @inbounds for i = 1:n-1
    κ₊ₕ[i] = κ[i] * κ[i+1] * (z[i] - z[i+1]) /
             (κ[i] * (z₊ₕ[i] - z[i+1]) + κ[i+1] * (z[i] - z₊ₕ[i])) # Eq. 5.16
  end

  a = zeros(n)
  b = zeros(n)
  c = zeros(n)
  d = zeros(n)

  # implicit
  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -κ₊ₕ[i] / dz₊ₕ[i]
      b[i] = cv[i] * dz[i] / dt - c[i] - df0
      d[i] = -κ₊ₕ[i] * (Tsoil_cur[i] - Tsoil_cur[i+1]) / dz₊ₕ[i] + f0
    elseif i < n
      a[i] = -κ₊ₕ[i-1] / dz₊ₕ[i-1]
      c[i] = -κ₊ₕ[i] / dz₊ₕ[i]
      b[i] = cv[i] * dz[i] / dt - a[i] - c[i]
      d[i] = κ₊ₕ[i-1] * (Tsoil_cur[i-1] - Tsoil_cur[i]) / dz₊ₕ[i-1] -
             κ₊ₕ[i] * (Tsoil_cur[i] - Tsoil_cur[i+1]) / dz₊ₕ[i]
    elseif i == n
      a[i] = -κ₊ₕ[i-1] / dz₊ₕ[i-1]
      c[i] = 0
      b[i] = cv[i] * dz[i] / dt - a[i]
      d[i] = κ₊ₕ[i-1] * (Tsoil_cur[i-1] - Tsoil_cur[i]) / dz₊ₕ[i-1]
    end
  end

  # --- Begin tridiagonal solution: forward sweep for layers N to 1
  # Bottom soil layer
  e = zeros(n)
  f = zeros(n)
  
  e[n] = a[n] / b[n]
  f[n] = d[n] / b[n]

  # Layers n-1 to 2
  @inbounds for i = n-1:-1:2
    den = b[i] - c[i] * e[i+1]
    e[i] = a[i] / den
    f[i] = (d[i] - c[i] * f[i+1]) / den
  end

  # Complete the tridiagonal solution to get the temperature of the top soil layer
  i = 1
  den = b[i] - c[i] * e[i+1]
  num = d[i] - c[i] * f[i+1]
  f[1] = num / den

  tsoi_test = Tsoil_cur[i] + f[1]

  # Potential melt rate based on temperature above freezing
  pot_snow_melt = max(0, (tsoi_test - tfrz) * den / λ_fus)
  max_snow_melt = snow_water / dt        # the amount of snow that is present
  snow_melt = min(max_snow_melt, pot_snow_melt) # Cannot melt more snow than present
  G_snow = snow_melt * λ_fus # Energy flux for snow melt

  # Update temperature
  Tsoil = zeros(n)
  Tsoil[1] = Tsoil_cur[1] + (num - G_snow) / den

  dtsoi = Tsoil[1] - Tsoil_cur[1]

  # Now complete the tridiagonal solution for layers 2 to N
  @inbounds for i = 2:n
    # dtsoi = f[i] - e[i] * Tsoil[i-1]
    dtsoi = f[i] - e[i] * dtsoi
    Tsoil[i] = Tsoil_cur[i] + dtsoi
  end
  # G_soil = f0([T_1]n) + df0 / dT * ([T_1]n + 1 - [T_1]n)
  G_soil = f0 + df0 * (Tsoil[1] - Tsoil_cur[1]) - G_snow # 一部分能量分配给融雪

  Tsoil, G_soil, G_snow
end


export soil_temperature_delta
