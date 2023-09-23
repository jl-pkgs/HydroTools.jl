"""
  PMLV2 (Penman–Monteith–Leuning Version 2) Evapotranspiration model

# Arguments

- `Prcp` : mm/d
- `Tavg` : degC
- `Rs`   : W m-2
- `Pa`   : kPa
- `U2`   : m/s
- `LAI`  : m2 m-2

# References
1. 
"""
function PMLV2(Prcp::T, Tavg, Rs, Rn, VPD, U2, LAI,
  Pa=atm, Ca=380.0; par) where {T<:Real}

  @unpack S_sls, kQ = par
  # D0 = 0.7
  # kQ = 0.6
  # kA = 0.7 # extinction coefficient
  PAR = 0.45 * Rs # W m-2, taken as 0.45 time of solar radiation
  PAR_mol = PAR * 4.57 # 1 W m-2 = 4.57 umol m-2 s-1

  λ, slope, gamma, Eeq = ET0_eq(Rn, Tair, Pa)
  ϵ = slope/gamma

  ### CARBON MODULE: PHOTOSYNTHESIS --------------------------------------------
  Am = par.Am_25 * fT2(Tavg)

  P1 = Am * par.Alpha * par.Thelta * PAR_mol
  P2 = Am * par.Alpha * PAR_mol
  P3 = Am * par.Thelta * Ca
  P4 = par.Alpha * par.Thelta * PAR_mol * Ca

  ## canopy conductance in (mol m-2 s-1)
  Ags = Ca * P1 / (P2 * kQ + P4 * kQ) * (
    kQ * LAI + log((P2 + P3 + P4) / (P2 + P3 * exp(kQ * LAI) + P4))) # umol m-1 s-1
  Ag = Ags  # gross assimilation rate in umol m-2 s-1
  Ag = Ag * f_VPD_Zhang2019 # * data$dhour_norm^2      # constrained by f_VPD;

  GPP = Ag * 86400 / 10^6 * 12 # [umol m-2 s-1] to [g C m-2 d-1]

  f_VPD_gc = 1 / (1 + VPD / D0) # Leuning f_vpd
  Gc = par.m * Ag / Ca * f_VPD_gc * 1.6 # g_water  = 1.6 * g_CO2 (mol m-2 s-1)
  
  ## Convert from mol m-2 s-1 to m s-1
  Gc = Gc * 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3)) # unit convert to m s-1
  Gc = max(Gc, 1e-6)
  
  ### WATER MODULE: ------------------------------------------------------------
  Ga = aerodynamic_conductance(U2, par.hc) # Leuning, 2008, Eq.13, doi:10.1029/2007WR006562

  # Transpiration from plant cause by radiation water transfer
  Tou = exp(-kA * LAI)
  LEcr = ϵ * Rn * (1 - Tou) / (ϵ + 1 + Ga / Gc) # W m-2

  # Transpiration from plant cause by aerodynamic water transfer
  ρ_a = cal_rho_a(Tair, q)
  LEca = (ρ_a * Cp * Ga * VPD / gamma) / (ϵ + 1 + Ga / Gc) # W m-2

  Ecr = W2mm(LEcr, λ) # [W m-2] change to [mm d-1]
  Eca = W2mm(LEca, λ) # [W m-2] change to [mm d-1]
  Ec = Ecr + Eca

  ## 5. Intercepted Evaporation (Ei)
  Ei = cal_Ei_Dijk2021(Prcp, LAI, par)
  Pi = Prcp - Ei

  Es_eq = Eeq * Tou # Soil evaporation at equilibrium, mm d-1
  GPP, Ec, Ei, Ecr, Eca, Ga, Gc, Pi, Es_eq
end





"""
    aerodynamic_conductance(U2, hc)

# Arguments
- `U2`: wind speed at 2m
- `hc`: canopy height

# Return
- `Ga`: aerodynamic conductance in m/s
"""
function aerodynamic_conductance(U2, hc)
  kmar = 0.40        # von Karman's constant 0.40
  d = 0.64 * hc
  zom = 0.13 * hc
  zoh = 0.10 * zom
  uz = cal_Uz(U2, Zob) # convert from u2 to uz
  Ga = uz * kmar^2 / (log((Zob - d) / zom) * log((Zob - d) / zoh)) # m s-1
  Ga
end


"""
# References
1. van Dijk, A.I.J.M, 2001, Eq2.
"""
function cal_Ei_Dijk2021(Prcp, LAI, par)
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


# 最大羧化能力温度调节函数
# V_m = Vm_25 * fT2
# `T` in -100:100, `fT2` always < 0.92
function fT2(Tavg::T)::T where {T<:Real}
  a = 0.031
  b = 0.115
  exp(a * (Tavg - 25.0)) / (1.0 + exp(b * (Tavg - 41.0))) # Gan2018, Eq. A5
end


# Piecewise function by Yongqiang and GanRong, 2019
function f_VPD_Zhang2019(VPD::T, par)::T where {T<:Real}
  if (VPD > par.VPDmax)
    T(0.0)
  elseif VPD < par.VPDmin
    T(1.0)
  else
    (VPDmax - VPD) / (VPDmax - VPDmin)
  end
end


export fT2, f_VPD_Zhang2019
