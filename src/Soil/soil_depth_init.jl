"""
    soil_depth_init(dz::AbstractVector)
    
Soil depth initialization

```julia
nsoil, z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
```
"""
function soil_depth_init(dz::AbstractVector)
  # Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
  # z_{i+1/2}
  nsoil = length(dz)

  z = zeros(nsoil)
  z₊ₕ = zeros(nsoil)
  dz₊ₕ = zeros(nsoil)

  z₊ₕ[1] = -dz[1]
  for i = 2:nsoil
    z₊ₕ[i] = z₊ₕ[i-1] - dz[i] # on the edge
  end

  # Soil depth (m) at center of layer i (negative distance from surface)
  z[1] = 0.5 * z₊ₕ[1]
  for i = 2:nsoil
    z[i] = 0.5 * (z₊ₕ[i-1] + z₊ₕ[i]) # on the center
  end

  # Thickness between between z(i) and z(i+1)
  for i = 1:nsoil-1
    dz₊ₕ[i] = z[i] - z[i+1]
  end
  dz₊ₕ[nsoil] = 0.5 * dz[nsoil]

  (; z, z₊ₕ, dz₊ₕ)
end
