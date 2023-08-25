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


MJ2W(x::Real) = x / 86400 * 1e6

MJ2mm(x::Real) = x / 2.45

W2MJ(x::Real) = x / 1e6 * 86400

function W2mm(x::Real; Tair=0)
  # Cp = 4.2 * 0.242  # specific heat at constant pressure, 1.013 [kJ kg⁻¹ °C⁻¹]
  lamada = 2500 - 2.2 * Tair
  return x / lamada * 86400 * 1e-3  # W M⁻² to mm
end


F2C(T_degF::Real) = (T_degF - 32.0) / (9.0 / 5.0)

C2F(T_degC::Real) = T_degC * (9.0 / 5.0) + 32

K2C(x::Real) = x - 273.15

C2K(x::Real) = x + 273.15


deg2rad(deg::Real) = deg / 180.0 * π

rad2deg(rad::Real) = rad / π * 180.0


export R2Q, Q2R, MJ2W, MJ2mm, W2MJ, W2mm, F2C, C2F, K2C, C2K, deg2rad, rad2deg
