# using Roots: find_zero
function find_zero_bisection(f, (a, b); tol=1e-6)
  mid = (a + b) / 2.0
  iters = 0
  while abs(f(mid)) > tol
    iters += 1
    if f(a) * f(mid) < 0
      b = mid
    else
      a = mid
    end
    mid = (a + b) / 2.0
  end
  return mid
end

"""
    adiabat_dry_T(P0::FT, Tk0::FT, P::FT, w=nothing)
    adiabat_dry_P(P0::FT, Tk0::FT, T::FT, w=nothing)

    LCL(T0::FT, Td::FT)
    LCL_bolton(T0::FT, Td::FT)

    theta(p::FT, tk::FT; p0=1000)
    theta_se(p, tk, w=0; p0=1000)
    
    theta_wet(P0, T0, Td)
    theta_wet_bolton(P0, T0, Td)
    
# Functions 

- `LCL`       : Bolton 1980, Eq. 15
- `LCL_bolton`: Bolton 1980, Eq. 22

# Arguments

- p: surface pressure, hPa
- Tk0, tk: surface temperature, K

- T0: surface temperature, degC
- P0: surface pressure, hPa
- Td: dew temperature, degC

# References
1. https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.equivalent_potential_temperature.html
2. https://github.com/wrf-model/WRF/blob/master/phys/module_diag_functions.F#L122
3. https://github.com/Unidata/MetPy/blob/e0e24d51702787943fc3c0481fa9a6632abe9d20/src/metpy/calc/thermo.py#L1512

# Examples
```julia
T1 = adiabat_dry_T(850., 20 + K0, 1000.) - K0
P1 = adiabat_dry_P(850., 20 +K0, T1+K0)

theta(850., 20.0 +K0)
theta_se(850., 20.0 +K0)

theta_wet(850., 20.0, 18.0)
theta_wet_bolton(850., 20.0, 18.0)
```
"""
function adiabat_dry_T(P0::FT, Tk0::FT, P::FT, w=nothing) where {FT<:Real}
  if w === nothing
    m = Rd / (Cp * 1e6) # 0.283
  else
    m = 0.2854 * (1 - 0.28 * w) # Bolton 1980
  end
  Tk0 * (P / P0)^m
end


function adiabat_dry_P(P0::FT, Tk0::FT, T::FT, w=nothing) where {FT<:Real}
  if w === nothing
    m = Rd / (Cp * 1e6) # 0.283
  else
    m = 0.2854 * (1 - 0.28 * w) # Bolton 1980
  end
  P0 * (T / Tk0)^(1 / m)
end


function theta(p::FT, tk::FT; p0=1000) where {FT<:Real}
  # R = 287.0    # J K-1 kg-1
  # cp = 1013 # J K-1 kg-1
  adiabat_dry_T(p, tk, FT(p0))
  # tk * (p0 / p)^(R / (Cp * 1e6))
end

function theta_se(p, tk, w=0; p0=1000)
  λ = cal_lambda(tk - K0) # MJ kg-1
  exp(λ * w / (Cp * tk)) * theta(p, tk; p0)
end

function LCL(T0::FT, Td::FT) where {FT<:Real}
  P0::FT = 1000
  ea = cal_es(Td) * 10
  w = ea2w(ea, P0) # g / g

  Tk0::FT = T0 + K0
  function goal(P)
    tc = adiabat_dry_T(P0, Tk0, P) - K0 # in Kdeg
    es = cal_es(tc) * 10
    ws = ea2w(es, P)
    w - ws
  end

  P_lcl = find_zero_bisection(goal, (20.0, P0)) # 1000 hPa to 20 hPa
  T_lcl = adiabat_dry_T(P0, Tk0, P_lcl) - K0
  T_lcl
end

function LCL_bolton(T0::FT, Td::FT) where {FT<:Real}
  RH = cal_es(Td) / cal_es(T0) * 100
  TK = T0 + K0
  term2 = 1 / (TK - 55) - log(RH / 100) / 2840
  1 / term2 + 55 - K0
end


function theta_wet(P0, T0, Td)
  T_lcl = LCL(T0, Td)
  K_lcl = T_lcl + K0

  Lv = cal_lambda(T_lcl)
  w = Tdew2w(Td, P0 / 10) # w守恒

  θ = adiabat_dry_T(P0, T0 + K0, 1000.0, w) # hPa
  θ_se = θ * exp(Lv * w / (Cp * K_lcl))
  (T_lcl, θ=θ - K0, θ_se=θ_se - K0)
end


function theta_wet_bolton(P0, T0, Td)
  td = Td + K0
  tk = T0 + K0

  ea = cal_es(Td) * 10
  r = ea2w(ea, P0)
  t_l = 56 + 1 / (1 / (td - 56) + log(tk / td) / 800)

  th_l = theta(P0 - ea, tk) * (tk / t_l)^(0.28 * r)
  θ_se = th_l * exp(r * (1 + 0.448 * r) * (3036 / t_l - 1.78))
  (; T_lcl=t_l - K0, θ_se=θ_se - K0)
end


export theta, theta_se, LCL, LCL_bolton,
  adiabat_dry_T, adiabat_dry_P,
  theta_wet, theta_wet_bolton
