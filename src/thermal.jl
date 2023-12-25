using Roots

"""
- tk: surface temperature, K
- p: surface pressure, hPa
"""
function theta(p, tk; p0=1000)
  R = 287.0    # J K-1 kg-1
  cp = 1004.0 # J K-1 kg-1
  # cp = 1013 # J K-1 kg-1
  # p0 = 1000 # hPa
  tk * (p0 / p)^(R / cp)
end

"""
- 
"""
function theta_se(p, tk, w=0; p0=1000)
  # tk = tc + K0 
  λ = 2.5e6 # J kg-1
  cp = 1004 # J K-1 kg-1
  exp(λ * w / (cp * tk)) * theta(tk, p; p0) - K0
end


"""
# Examples
```julia
T1 = adiabat_dry_T(850., 20 + K0, 1000.) - K0
P1 = adiabat_dry_P(850., 20 +K0, T1+K0)
```
"""
function adiabat_dry_T(P₀::FT, Tk₀::FT, P::FT, w=nothing) where {FT<:Real}
  if w === nothing
    m = Rd / (Cp * 1e6) # 0.283
  else
    m = 0.2854 * (1 - 0.28 * w) # Bolton 1980
  end
  Tk₀ * (P / P₀)^m
end

function adiabat_dry_P(P₀::FT, Tk₀::FT, T::FT, w=nothing) where {FT<:Real}
  if w === nothing
    m = Rd / (Cp * 1e6) # 0.283
  else
    m = 0.2854 * (1 - 0.28 * w) # Bolton 1980
  end
  P₀ * (T / Tk₀)^(1 / m)
end


function LCL(T0::FT, Td::FT) where {FT<:Real}
  P0::FT = 1000
  ea = cal_es(Td) * 10
  w = ea2w(ea, P0) # g / g

  Tk₀::FT = T0 + K0
  function goal(P)
    tc = adiabat_dry_T(P0, Tk₀, P) - K0 # in Kdeg
    es = cal_es(tc) * 10
    ws = ea2w(es, P)
    w - ws
  end

  P_lcl = find_zero(goal, (20, P0)) # 1000 hPa to 20 hPa
  T_lcl = adiabat_dry_T(P0, Tk₀, P_lcl) - K0
  T_lcl
end

"""
- `LCL`       : Bolton 1980, Eq. 15
- `LCL_bolton`: Bolton 1980, Eq. 22
"""
function LCL_bolton(T0::FT, Td::FT) where {FT<:Real}
  RH = cal_es(Td) / cal_es(T0) * 100
  TK = T0 + K0
  term2 = 1 / (TK - 55) - log(RH / 100) / 2840
  1 / term2 + 55 - K0
end


"""
  theta_wet(P0, T0, Td)

# References
1. https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.equivalent_potential_temperature.html
2. https://github.com/wrf-model/WRF/blob/master/phys/module_diag_functions.F#L122
3. https://github.com/Unidata/MetPy/blob/e0e24d51702787943fc3c0481fa9a6632abe9d20/src/metpy/calc/thermo.py#L1512
"""
function theta_wet(P0, T0, Td)
  T_lcl = LCL(T0, Td)
  K_lcl = T_lcl + K0

  Lv = cal_lambda(T_lcl)
  w = Tdew2w(Td, P0 / 10) # w守恒

  Θ = adiabat_dry_T(P0, T0 + K0, 1000.0, w) # hPa
  Θ_se = Θ * exp(Lv * w / (Cp * K_lcl))
  (T_lcl, Θ=Θ - K0, Θ_se=Θ_se - K0)
end


function theta_wet_bolton(P0, T0, Td)
  ea = cal_es(Td) * 10

  td = Td + K0
  tk = T0 + K0

  r = ea2w(ea, P0)

  t_l = 56 + 1 / (1 / (td - 56) + log(tk / td) / 800)
  th_l = theta(P0 - ea, tk) * (tk / t_l)^(0.28 * r)
  Θ_se = th_l * exp(r * (1 + 0.448 * r) * (3036 / t_l - 1.78))
  (; T_lcl=t_l - K0, Θ_se=Θ_se - K0)
end


export theta, theta_se, LCL, LCL_bolton,
  adiabat_dry_T, adiabat_dry_P,
  theta_wet, theta_wet_bolton
