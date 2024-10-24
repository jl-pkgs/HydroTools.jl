export soil_moisture!

# TODO: 
# - ψ0没有参与更新，要如何解决？
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

## Arguments
- `θ`     : [m3 m-3]
- `ψ`     : [cm]，为负，约干约负
- `ψ0`    : 已知Q0，可以求得ψ0, Eq. 8.25
- `dz`    : [cm]
- `dt`    : [s]
- `param` : 参数

- `sink`  : 蒸发项，[cm s-1]
- `Q0`    : 下渗，[cm s-1]

# TODO: 
1. 加入sink相：体现蒸发`E`和入渗`I`的影响（I的影响，主要发生在前两层）
2. 更新每一层的土壤深度

# Example
```julia
soil_moisture!(θ, ψ, ψ0, dz, dt, param)
```
"""
function soil_moisture!(θ, ψ, ψ0, dz, dt, param; sink=nothing, Q0=nothing, fun=van_Genuchten)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  n = length(dz)

  # a deep copy
  θ_n = deepcopy(θ)
  ψ_n = deepcopy(ψ)
  isnothing(sink) && (sink = zeros(n))

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

  K0₊ₕ = K[1]
  dz0₊ₕ = 0.5 * dz[1]
  !isnothing(Q0) && (ψ0 = -(Q0 / K0₊ₕ + 1) * dz0₊ₕ + ψ[1]) # Eq. 8.25

  a = zeros(n)
  c = zeros(n)
  b = zeros(n)
  d = zeros(n)

  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / dz₊ₕ[i]
      b[i] = Cap[i] * dz[i] / (0.5 * dt) + K0₊ₕ / dz0₊ₕ - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
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
    d[i] -= sink[i]
  end
  ψ_pred .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1/2 time

  ## update: θ, K and Cap
  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ_pred[i]; param)
  end
  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
  K0₊ₕ = K[1] # 可以按照同样的方法，设置
  !isnothing(Q0) && (ψ0 = -(Q0 / K0₊ₕ + 1) * dz0₊ₕ + ψ[1]) # Eq. 8.25

  ## second round
  # Terms for tridiagonal matrix
  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * dz₊ₕ[i])
      b[i] = Cap[i] * dz[i] / dt - c[i] + K0₊ₕ / (2 * dz0₊ₕ)
      d[i] = Cap[i] * dz[i] / dt * ψ[i] +
             K0₊ₕ / (2 * dz0₊ₕ) * (2ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
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
    d[i] -= sink[i]
  end
  ψ .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1

  # --- Check water balance
  Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_n[1]) + (ψ0 - ψ[1])) - K0₊ₕ
  QN = -K[n]

  dθ = 0
  for i = 1:n
    dθ += (θ[i] - θ_n[i]) * dz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
