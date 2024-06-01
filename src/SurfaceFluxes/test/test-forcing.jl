cal_ρ(RH2ea(100), 0)

# z = 0.0
Ta = 20.0
ea = RH2ea(60.) * 1e3
Pa = 95.0 * 1e3 # hPa
z = 2.0

# RH = ea / (cal_es(Ta)*10) * 100
# cpₐ, T_pot, Tv, mm_air, ρ_air, ρ_mol = cal_Cp(Ta, ea, Pa, z)
# ρₐ = cal_ρ(Ta, ea, Pa)
init_forcing(Ta, ea, Pa, z)
