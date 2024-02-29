export soil_moisture!

# TODO: 
# - ψ0没有参与更新，要如何解决？
# - 下渗，Q0
# - 添加一个结构体，保存中间变量，节省内存


"""
    soil_moisture!(θ, ψ, ψ0, dz, dt, param; fun=van_Genuchten)

为了解决相互依赖的问题，这里求解转了两次。
```julia
# 计算n+1/2时刻的ψ_pred，explicit方案：空间上采用n-1时刻的差分
θ, K, Cap -> ψ_pred

# Crank-Nicolson方案：空间上采用n+1, n-1时刻
ψ_pred -> θ, K, Cap
θ, K, Cap -> ψ
```

# Example
```julia
soil_moisture!(Θ, ψ, ψ0, dz, dt, param)
```
"""
function soil_moisture!(θ, ψ, ψ0, dz, dt, param; fun=van_Genuchten)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  n = length(dz)

  # a deep copy
  θ_n = deepcopy(θ)
  ψ_n = deepcopy(ψ)

  # θ = zeros(n)
  K = zeros(n)
  Cap = zeros(n)
  ψ_pred = zeros(n)
  K₊ₕ = zeros(n)

  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  end
  
  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end

  K0_₊ₕ = K[1]
  dz0_₊ₕ = 0.5 * dz[1]

  a = zeros(n)
  c = zeros(n)
  b = zeros(n)
  d = zeros(n)

  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / dz₊ₕ[i]
      b[i] = Cap[i] * dz[i] / (0.5 * dt) + K0_₊ₕ / dz0_₊ₕ - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K0_₊ₕ / dz0_₊ₕ * ψ0 + K0_₊ₕ - K₊ₕ[i]
    elseif i < n
      a[i] = -K₊ₕ[i-1] / dz₊ₕ[i-1]
      c[i] = -K₊ₕ[i] / dz₊ₕ[i]
      b[i] = Cap[i] * dz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == n
      a[i] = -K₊ₕ[n-1] / dz₊ₕ[n-1]
      c[i] = 0
      b[i] = Cap[i] * dz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[n-1] - K[i]
    end
  end
  ψ_pred .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1/2 time

  ## update: θ, K and Cap
  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ_pred[i]; param)
  end
  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
  K0_₊ₕ = K[1]
  
  ## second round
  # Terms for tridiagonal matrix
  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * dz₊ₕ[i])
      b[i] = Cap[i] * dz[i] / dt + K0_₊ₕ / (2 * dz0_₊ₕ) - c[i]
      d[i] = Cap[i] * dz[i] / dt * ψ[i] + K0_₊ₕ / (2 * dz0_₊ₕ) * ψ0 +
             K0_₊ₕ / (2 * dz0_₊ₕ) * (ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0_₊ₕ - K₊ₕ[i]
    elseif i < n
      a[i] = -K₊ₕ[i-1] / (2 * dz₊ₕ[i-1])
      c[i] = -K₊ₕ[i] / (2 * dz₊ₕ[i])
      b[i] = Cap[i] * dz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    else
      i == n
      a[i] = -K₊ₕ[i-1] / (2 * dz₊ₕ[i-1])
      c[i] = 0
      b[i] = Cap[i] * dz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
  end
  ψ .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1

  # --- Check water balance
  Q0 = -K0_₊ₕ / (2 * dz0_₊ₕ) * ((ψ0 - ψ_n[1]) + (ψ0 - ψ[1])) - K0_₊ₕ
  QN = -K[n]

  dθ = 0
  for i = 1:n
    dθ += (θ[i] - θ_n[i]) * dz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
