export Norman_Longwave;

"""
    Norman_Longwave(ϵ=0.98, ϵ_g=1.0, T_veg=25, T_g=20, L_sky=400)

    Longwave radiation transfer through canopy using Norman (1979)

# Arguments
- `ϵ`     : Leaf emissivity
- `ϵ_g`   : Ground (soil) emissivity
- `LAI`   : Leaf area index (m2/m2)
- `L_sky` : Atmospheric longwave radiation (W/m2)

# Examples
```julia
U, L_up, L_dn = Norman_Longwave()
```
"""
function Norman_Longwave(ϵ=0.98, ϵ_g=1.0, T_veg=25, T_g=20, L_sky=400; check_error=true)
  # 这里应该定义一个结构体，存放参数
  nveg = 49                    # Number of leaf layers (each with lai = dlai)
  nsoi = 1                     # First canopy layer is soil
  nbot = nsoi + 1              # Index for bottom leaf layer
  ntop = nbot + nveg - 1       # Index for top leaf layer, 50

  tleaf = ones(ntop) .* T_veg  # Leaf temperature (degC)
  τd = ones(ntop) .* 0.915     # transmittance of diffuse radiation through a single leaf layer
  # dlai = ones(ntop) .* 0.1   # Layer leaf area index (m2/m2)

  ω = 1 - ϵ # Leaf scattering coefficient
  ρ = ω     # Leaf reflectance
  τ = 0     # Leaf transmittance
  
  params = (; tleaf, τd, τ, ρ, ϵ, ϵ_g)
  # @show params

  ## 1. solve longwave radiation -----------------------------------------------
  nlayers = ntop
  atri = zeros(nlayers * 2)
  btri = zeros(nlayers * 2)
  ctri = zeros(nlayers * 2)
  dtri = zeros(nlayers * 2)

  # Emitted longwave radiation from leaves (W/m2)
  ir_source = blackbody.(tleaf, ϵ) .* (1 .- τd) # outward Longwave

  # global variable: atri, btri, ctri, dtri, ir_source, τd
  function update_ef(i, m; direction="up")
    refld = (1 - τd[i]) * ρ
    trand = (1 - τd[i]) * τ + τd[i]

    if direction == "up"
      fiv = refld - trand * trand / refld
      eiv = trand / refld

      atri[m] = -eiv
      btri[m] = 1
      ctri[m] = -fiv
      dtri[m] = (1 - eiv) * ir_source[i]
    elseif direction == "down"
      aiv = refld - trand * trand / refld
      biv = trand / refld
      atri[m] = -aiv
      btri[m] = 1
      ctri[m] = -biv
      dtri[m] = (1 - biv) * ir_source[i]
    end
  end

  # Soil: upward flux
  iv = nsoi
  m = 1
  atri[m] = 0
  btri[m] = 1
  ctri[m] = -(1 - ϵ_g)
  dtri[m] = blackbody(T_g, ϵ_g)

  # Soil: downward flux
  update_ef(iv + 1, 2; direction="down")

  # Leaf layers, excluding top layer
  for iv = nbot:ntop-1
    i1 = (iv - 1) * 2 + 1
    i2 = (iv - 1) * 2 + 2
    update_ef(iv, i1; direction="up")       # upward flux
    update_ef(iv + 1, i2; direction="down") # downward flux
  end

  # Top canopy laye: upward flux
  iv = ntop
  i1 = (iv - 1) * 2 + 1
  i2 = (iv - 1) * 2 + 2
  update_ef(iv, i1)
  # Top canopy layer: downward flux
  m = i2
  atri[m] = 0
  btri[m] = 1
  ctri[m] = 0
  dtri[m] = L_sky

  # Solve tridiagonal equations for upward and downward fluxes
  U = tridiagonal_solver(atri, btri, ctri, dtri, m) # Eq. 124
  L_up, L_dn = U2longwave(U; nsoi, nbot) # Soil fluxes

  L_g = blackbody(T_g, ϵ_g)
  check_error && Norman_checkError(L_up, L_dn, L_sky, L_g, params; nsoi, ntop)

  # Rln for each layer
  Rln = zeros(nlayers)
  for i = nsoi+1:ntop
    L = blackbody(tleaf[i], ϵ)
    Rln[i] = (ϵ * (L_dn[i] + L_up[i-1]) - 2 * L) * (1 - τd[i]) # Eq. 14.129
  end
  Rln[1] = L_dn[1] - L_up[1]
  e = sum(Rln) - (L_sky - L_up[end]) # Eq. 14.333
  @assert abs(e) <= 1e-8

  Rln_soil = Rln[1]
  Rln_veg = sum(Rln[2:end])

  (; L_up, L_dn, Rln, Rln_soil, Rln_veg)
end


# L_up, L_dn = U2longwave(U; nsoi)
function U2longwave(U; nsoi=1, nbot)
  nlayers = floor(Int, length(U) / 2)
  ntop = nlayers
  L_up = zeros(nlayers)
  L_dn = zeros(nlayers)

  # Soil fluxes
  i = nsoi
  L_up[i] = U[1]
  L_dn[i] = U[2]
  # Leaf layer fluxes
  for i = nbot:ntop
    i1 = (i - 1) * 2 + 1
    i2 = (i - 1) * 2 + 2
    L_up[i] = U[i1]
    L_dn[i] = U[i2]
  end
  L_up, L_dn
end

# 检验矩阵求解的正确性
# Error check: compare tridiagonal solution with actual equations
function Norman_checkError(L_up, L_dn, L_sky, L_g, params; nsoi, ntop)
  @unpack tleaf, τd, τ, ρ, ϵ, ϵ_g = params

  nlayers = length(L_up)
  L_dn_eq = zeros(nlayers)
  L_up_eq = zeros(nlayers)

  i = ntop
  L_dn_eq[i] = L_sky
  for i = ntop-1:-1:nsoi
    L_dn_eq[i] = L_dn[i+1] * (τd[i+1] + (1 - τd[i+1]) * τ) +
                 L_up[i] * ((1 - τd[i+1]) * ρ) +
                 blackbody(tleaf[i+1], ϵ) * (1 - τd[i+1]) # Bonan2019, Eq. 14.122
  end

  i = nsoi
  L_up_eq[i] = (1 - ϵ_g) * L_dn[i] + L_g
  for i = nsoi:ntop-1
    L_up_eq[i+1] = L_up[i] * (τd[i+1] + (1 - τd[i+1]) * τ) +
                   L_dn[i+1] * ((1 - τd[i+1]) * ρ) +
                   blackbody(tleaf[i+1], ϵ) * (1 - τd[i+1]) # Bonan2019, Eq. 14.123
  end

  for i = nsoi:ntop
    err = L_dn_eq[i] - L_dn[i]
    if (abs(err) > 1e-10)
      @printf("err = %15.5f\n", err)
      @printf("tridiag = %15.5f\n", L_dn[i])
      @printf("eq = %15.5f\n", L_dn_eq[i])
      error("Norman radiation: downward error")
    end
    err = L_up_eq[i] - L_up[i]
    if (abs(err) > 1e-10)
      @printf("err = %15.5f\n", err)
      @printf("tridiag = %15.5f\n", L_up[i])
      @printf("eq = %15.5f\n", L_up_eq[i])
      error("Norman radiation: upward error")
    end
  end
end
