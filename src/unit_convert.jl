"""
    R2Q(R [, area = dt * 3.6, dt = 24])
    Q2R(Q [, area = dt * 3.6, dt = 24])

Convert runoff from mm per dt to flow rate in m³/s.

## Arguments
- `R::Real`: Runoff (mm per dt).
- `Q::Real`: Flow rate in m³/s.
- `area::Real`: Basin area in km². Default is `dt * 3.6`.
- `dt::Real`: Time duration in hours. Default is 24.

## Returns
Flow rate in m³/s.

"""
R2Q(R::Real, area::Real; dt=24) = R * area / (dt * 3.6)

Q2R(Q::Real, area::Real; dt=24) = Q * dt * 3.6 / area


MJ2W(x::Real) = x / 0.0864 # x / 86400 * 1e6

MJ2mm(x::Real) = x / 2.45

W2MJ(x::Real) = x * 0.0864

function W2mm(x::Real, Tair::Real)
  lamada = 2500 - 2.2 * Tair
  return x / lamada * 86400 * 1e-3  # W M⁻² to mm
end

# lambda: [MJ kg-1]
W2mm(Ra; lambda) = Ra * 86400 / 1e6 / lambda


"""
    mol2m(Tavg, Pa=atm)
    mol2m_rong2018(Tavg, Pa=atm)

Convert from mol m-2 s-1 to m s-1, g = g_m * mol2m(Tavg)

# Reference
1. Monteith, 2013, Principles of Environmental Physics, Eq. 3.14
"""
mol2m(Tavg, Pa=atm) = R * (Tavg + K0) / Pa / 1000

mol2m_rong2018(Tavg, Pa=atm) = 1e-2 / (0.446 * (273 / (273 + Tavg)) * (Pa / 101.3))


function F2C(T_degF::FT)::FT where {FT<:Real} 
  (T_degF - 32) / (9.0 / 5.0)
end

function C2F(T_degC::FT)::FT where {FT<:Real}
  T_degC * (9.0 / 5.0) + 32
end

K2C(x::Real) = x - 273.15

C2K(x::Real) = x + 273.15


deg2rad(deg::Real) = deg / 180.0 * π

rad2deg(rad::Real) = rad / π * 180.0


export R2Q, Q2R, MJ2W, MJ2mm, W2MJ, W2mm, F2C, C2F, K2C, C2K
# deg2rad, rad2deg
