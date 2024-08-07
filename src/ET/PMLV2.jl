# 相同植被类型多个站点一起的计算
function PMLV2_sites(df::AbstractDataFrame; par::Param_PMLV2=Param_PMLV2(), kw...)
  sites = df.site
  grps = unique(sites)

  # ! 这里有优化空间，但程序会写的很复杂
  res = []
  for grp in grps
    inds = sites .== grp
    d = df[inds, :]
    r = PMLV2(d; par, kw...)
    push!(res, r)
  end
  vcat(res...)
end


"""
    PMLV2(Prcp, Tavg, Rs, Rn, VPD, U2, LAI, Pa, Ca; par=param0, frame=3)

# Notes
一个站点的计算。注意，不同植被类型，参数不同。

# Arguments
- `frame`: in 8-days
"""
function PMLV2(Prcp::V, Tavg::V, Rs::V, Rn::V,
  VPD::V, U2::V, LAI::V,
  Pa::V, 
  Ca::Union{T,V}=T(380.0);
  par::Param_PMLV2=Param_PMLV2(), frame=3,
  res::Union{Nothing,output_PML}=nothing) where {T<:Real, V<:AbstractVector{T}}

  n = length(Prcp)
  fields = fieldnames(interm_PML)[2:end-2] # skip ET, fval_soil and Es
  r = interm_PML{T}()
  res === nothing && (res = output_PML{T}(; n))

  isvec_Ca = isa(Ca, AbstractVector)
  @inbounds for t = 1:n
    _Ca = isvec_Ca ? Ca[t] : Ca

    PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t], 
      _Ca; par, r)
    res[t, fields] = r
  end

  res.fval_soil = movmean2(Prcp, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.fval_soil, 0, 1)

  res.Es .= res.fval_soil .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end


"""
# TODO
> 

# Arguments
- `kw`: named keyword arguments
  + `r`: `interm_PML`
"""
function PMLV2(d::AbstractDataFrame; par::Param_PMLV2=Param_PMLV2(), kw...)
  PMLV2(d.Prcp, d.Tavg, d.Rs, d.Rn,
    d.VPD, d.U2, d.LAI,
    d.Pa, d.Ca; par, kw...) |> to_df
end


"""
  PMLV2 (Penman–Monteith–Leuning Version 2) Evapotranspiration model

# Arguments

- `Prcp` : mm/d
- `Tavg` : degC
- `Rs`   : W m-2
- `Rn`   : W m-2
- `VPD`  : W m-2
- `U2`   : m/s
- `LAI`  : m2 m-2
- `Pa`   : kPa
- `Ca`   : ppm

# Examples
```julia
```

# References
1. Gan Rong, 2018, Ecohydrology
2. Zhang Yongqiang, 2019, RSE
3. Kong Dongdong, 2019, ISPRS
"""
function PMLV2(Prcp::T, Tavg::T, Rs::T, Rn::T, VPD::T, U2::T, LAI::T,
  Pa=atm, Ca=380.0; 
  par::Param_PMLV2=Param_PMLV2(),
  r::Union{Nothing,interm_PML}=nothing) where {T<:Real}
  r === nothing && (r = interm_PML{T}())

  # D0 = 0.7
  # kQ = 0.6 # extinction coefficients for visible radiation
  # kA = 0.7 # extinction coefficient for available energy
  λ, Δ, γ, r.Eeq = ET0_eq(Rn, Tavg, Pa)
  ϵ = Δ / γ

  ### CARBON MODULE: PHOTOSYNTHESIS --------------------------------------------
  r.GPP, r.Gc_w = photosynthesis(Tavg, Rs, VPD, LAI, Pa, Ca; par)

  ### WATER MODULE: ------------------------------------------------------------
  r.Ga = aerodynamic_conductance(U2, par.hc) # Leuning, 2008, Eq.13, doi:10.1029/2007WR006562

  # Transpiration from plant cause by radiation water transfer
  Tou = exp(-par.kA * LAI)
  LEcr = ϵ * Rn * (1 - Tou) / (ϵ + 1 + r.Ga / r.Gc_w) # W m-2

  # Transpiration from plant cause by aerodynamic water transfer
  ρ_a = cal_rho_a(Tavg, Pa)
  # ρ_a = cal_rho_a(Tavg, q, Pa)
  LEca = (ρ_a * Cp * 1e6 * r.Ga * VPD / γ) / (ϵ + 1 + r.Ga / r.Gc_w) # W m-2, `Cp*1e6`: [J kg-1 °C-1]

  r.Ecr = W2mm(LEcr; lambda=λ) # [W m-2] change to [mm d-1]
  r.Eca = W2mm(LEca; lambda=λ) # [W m-2] change to [mm d-1]
  r.Ec = r.Ecr + r.Eca

  ## 5. Intercepted Evaporation (Ei)
  r.Ei = cal_Ei_Dijk2021(Prcp, LAI, par)
  r.Pi = Prcp - r.Ei

  ## TODO: 补充冰面蒸发的计算
  Evp::T = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ
  r.ET_water = r.Eeq + Evp

  r.Es_eq = r.Eeq * Tou # Soil evaporation at equilibrium, mm d-1
  r
  # GPP, Ec, Ecr, Eca, Ei, Pi, Es_eq, Eeq, ET_water, Ga, Gc_w
end



"""
    photosynthesis(Tavg::T, Rs::T, VPD::T, LAI::T, Ca=380.0; par)

# Example
```julia
# GPP, Gc_w = photosynthesis(Tavg, Rs, VPD, LAI, Ca; par)
```
"""
function photosynthesis(Tavg::T, Rs::T, VPD::T, LAI::T,
  Pa=atm, Ca=380.0; par::Param_PMLV2) where {T<:Real}

  kQ = par.kQ # light extinction coefficient

  PAR = 0.45 * Rs # W m-2, taken as 0.45 time of solar radiation
  PAR_mol = PAR * 4.57 # 1 W m-2 = 4.57 umol m-2 s-1

  Vm = par.Am_25 * T_adjust_Vm25(Tavg)
  Am = Vm # 认为最大光合速率 = 最大羧化能力

  P1 = Am * par.Alpha * par.Thelta * PAR_mol
  P2 = Am * par.Alpha * PAR_mol
  P3 = Am * par.Thelta * Ca
  P4 = par.Alpha * par.Thelta * PAR_mol * Ca

  ## canopy conductance in (mol m-2 s-1)
  Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
    kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4))) # umol m-1 s-1
  Ag = Ags  # gross assimilation rate in umol m-2 s-1
  Ag = Ag * f_VPD_Zhang2019(VPD, par) # * data$dhour_norm^2  # constrained by f_VPD;

  GPP = Ag * 86400 / 10^6 * 12 # [umol m-2 s-1] to [g C m-2 d-1]

  f_VPD_gc = 1.0 / (1.0 + VPD / par.D0) # Leuning f_vpd
  Gc = par.m * Ag / Ca * f_VPD_gc # canopy conductance for carbon

  ## Convert from mol m-2 s-1 to m s-1
  Gc = Gc * 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3)) # Gc = Gc * mol2m(Tavg, Pa)
  Gc = max(Gc, 1e-6)

  Gc_w = Gc * 1.6 # g_water  = 1.6 * g_CO2 (mol m-2 s-1), canopy conductance for water
  GPP, Gc_w
end


# 最大羧化能力温度调节函数
# V_m = Vm_25 * T_adjust_Vm25
# `T` in -100:100, `T_adjust_Vm25` always < 0.92
function T_adjust_Vm25(Tavg::T)::T where {T<:Real}
  a = 0.031
  b = 0.115
  exp(a * (Tavg - 25.0)) / (1.0 + exp(b * (Tavg - 41.0))) # Gan2018, Eq. A5
end


# Piecewise function by Yongqiang and GanRong, 2019
function f_VPD_Zhang2019(VPD::T, par::Param_PMLV2)::T where {T<:Real}
  if (VPD > par.VPDmax)
    T(0.0)
  elseif VPD < par.VPDmin
    T(1.0)
  else
    (par.VPDmax - VPD) / (par.VPDmax - par.VPDmin)
  end
end


"""
# References
1. van Dijk, A.I.J.M, 2001, Eq2.
"""
function cal_Ei_Dijk2021(Prcp::T, LAI::T, par::Param_PMLV2) where {T<:Real}
  # two params in Ei
  # @unpack S_sls, fER0 = par
  LAIref = 5
  # van Dijk, A.I.J.M, 2001, Eq2.
  fveg = 1 - exp(-LAI / LAIref)  # Canopy cover fraction, Eq.1
  Sveg = par.S_sls * LAI # Specific leaf storage, Eq.2
  fER = par.fER0 * fveg # the value of 0.50 based on optimisation at Australian catchments
  Pwet = -log(1 - par.fER0) / par.fER0 * Sveg / fveg # -log(1 - fER /fveg),

  # Pwet[is.na(Pwet)] = 0; check, negative value, log will give error
  Ei = Prcp < Pwet ? fveg * Prcp : (fveg * Pwet + fER * (Prcp - Pwet))
  Ei
end


export PMLV2_sites
export PMLV2, T_adjust_Vm25, f_VPD_Zhang2019
