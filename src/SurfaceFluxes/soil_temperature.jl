"""
    soil_temperature(dz::AbstractVector, dt, 
        κ::AbstractVector, cv::AbstractVector, 
        Tsoil_cur::AbstractVector, Tsurf_next::Real;
        solution="implicit", method="apparent-heat-capacity")

> 这里采用的是方案1的边界条件, sum(F) = G, G = - k1 * (T1 - T0) / dz1

# Arguments: 

- `κ`: thermal conductivity (W/m/K)
- `cv`: volumetric heat capacity (J/m3/K)

- `Tsurf_next`: Tsurf_next_next, T0_{n+1}
- `solution`:
  + `implicit`:
  + `crank-nicolson`:

# Examples

```julia
Tsoil_next, G = soil_temperature(dz, dt, κ, cv, Tsoil_cur, Tsurf_next)
```
"""
function soil_temperature(dz::AbstractVector, dt, 
  κ::AbstractVector, cv::AbstractVector, 
  Tsoil_cur::AbstractVector, Tsurf_next::Real;
  solution="implicit", method="apparent-heat-capacity")

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

  ## Set up tridiagonal matrix
  @inbounds if solution == "implicit"
    for i = 1:n
      m = cv[i] * dz[i] / dt
      if i == 1
        a[i] = 0
        c[i] = -κ₊ₕ[i] / dz₊ₕ[i]
        b[i] = m - c[i] + κ[i] / (0 - z[i])  # κ_(1/2) / dz_(1/2) = κ₁ / (0 - z₁)
        d[i] = m * Tsoil_cur[i] + κ[i] / (0 - z[i]) * Tsurf_next

      elseif i < n
        a[i] = -κ₊ₕ[i-1] / dz₊ₕ[i-1]
        c[i] = -κ₊ₕ[i] / dz₊ₕ[i]
        b[i] = m - a[i] - c[i]
        d[i] = m * Tsoil_cur[i]

      elseif i == n
        a[i] = -κ₊ₕ[i-1] / dz₊ₕ[i-1]
        c[i] = 0
        b[i] = m - a[i]
        d[i] = m * Tsoil_cur[i]
      end
    end

  elseif solution == "crank-nicolson"
    # --- Heat flux at time n (W/m2) of each layer
    f = zeros(n)
    for i = 1:n-1
      f[i] = -κ₊ₕ[i] / dz₊ₕ[i] * (Tsoil_cur[i] - Tsoil_cur[i+1]) # Eq. 5.15
    end

    for i = 1:n
      m = cv[i] * dz[i] / dt
      if i == 1
        a[i] = 0
        c[i] = -0.5 * κ₊ₕ[i] / dz₊ₕ[i]
        b[i] = m - c[i] + κ[i] / (0 - z[i])
        d[i] = m * Tsoil_cur[i] + κ[i] / (0 - z[i]) * Tsurf_next + 0.5 * f[i]

      elseif i < n
        a[i] = -0.5 * κ₊ₕ[i-1] / dz₊ₕ[i-1]
        c[i] = -0.5 * κ₊ₕ[i] / dz₊ₕ[i]
        b[i] = m - a[i] - c[i]
        d[i] = m * Tsoil_cur[i] + 0.5 * (f[i] - f[i-1])

      elseif i == n
        a[i] = -0.5 * κ₊ₕ[i-1] / dz₊ₕ[i-1]
        c[i] = 0
        b[i] = m - a[i]
        d[i] = m * Tsoil_cur[i] - 0.5 * f[i-1]
      end
    end
  end

  Tsoil_next = tridiagonal_solver(a, b, c, d) # the updated Tsoil

  # --- Derive energy flux into soil (W/m2)
  G = κ[1] * (Tsurf_next - Tsoil_next[1]) / (0 - z[1]) # gound heat flux, G

  ## --- Check for energy conservation
  if method == "apparent-heat-capacity"
    LE_f = 0.0
  elseif method == "excess-heat"
    # phase_change(physcon, soilvar, dt);
  end

  edif = 0.0 # Sum change in energy (W/m2)
  @inbounds for i = 1:n
    edif += cv[i] * dz[i] * (Tsoil_next[i] - Tsoil_cur[i]) / dt
  end

  # Error check
  err = edif - G - LE_f
  abs(err) > 1e-03 && error("Soil temp erature energy conservation error")

  Tsoil_next, G
end
