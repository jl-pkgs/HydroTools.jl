using HydroTools
using Test

@test cal_ρ(RH2ea(100), 0) ≈ 1.276932763678108

# z = 0.0
Ta = 20.0
ea = RH2ea(60.) * 1e3
Pa = 95.0 * 1e3 # kPa
z = 2.0

init_forcing(Ta, ea, Pa, z)

# RH = ea / (cal_es(Ta)*10) * 100
# cpₐ, T_pot, Tv, mm_air, ρ_air, ρ_mol = cal_Cp(Ta, ea, Pa, z)
# ρₐ = cal_ρ(Ta, ea, Pa)


Mw = 18.01528
Md = 28.96340

M_dry = 18.01528 * 1e-3 # [kg mol-1]
M_h2o = 28.96340 * 1e-3

R = 8.3144621 # J/(mol K)

# Rw = R / M_h2o # ≈ 461.5 [J kg-1 K-1]
# Rd = R / M_dry # ≈ 287.1 [J kg-1 K-1]

Rw = R / Mw * 1000
Rd = R / Md * 1000


Mw = R / M_h2o
Md = R / M_dry
