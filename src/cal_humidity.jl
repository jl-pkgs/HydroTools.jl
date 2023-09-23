"""
    cal_humidity

# Arguments
- `w`    : mix ratio, m_w / m_d
- `q`    : specific humidity in kg/kg or g/g
- `Pa`   : surface air pressure
- `ea`   : actual vapor pressure (kPa)
- `RH`   : relative humidity, in %
- `Tair` : air temperature, in degC
- `Tdew` : dew temperature (in degC)

"""
function cal_es(Tair::T)::T where {T<:Real}
  0.6108 * exp((17.27 * Tair) / (Tair + 237.3)) # FAO98 equation
  # 0.61094 * exp((17.625 * Tair) / (Tair + 243.04))
end

function cal_es(Tmin::T, Tmax::T)::T where {T<:Real}
  (cal_es(Tmin) + cal_es(Tmax)) / T(2.0)
end

function Tdew2RH(Tdew::T, Tair::T)::T where {T<:Real}
  ea::T = cal_es(Tdew)
  es::T = cal_es(Tair)
  ea / es * 100
end

function Tdew2VPD(Tdew::T, Tair::T)::T where {T<:Real}
  ea::T = cal_es(Tdew)
  es::T = cal_es(Tair)
  max(es - ea, 0.0)
end

function ea2q(ea, Pa=atm)
  ϵ * ea / (Pa - (1 - ϵ) * ea)
end

function q2ea(q, Pa=atm)
  q * Pa/(ϵ + (1 - ϵ) * q)
end

function q2RH(q, Tair, Pa=atm)
  ea = q2ea(q, Pa)
  es = cal_es(Tair)
  ea / es * 100
end

function q2VPD(q, Tmin, Tmax, Pa=atm)
  ea = q2ea(q, Pa)
  es = cal_es(Tmin, Tmax)
  max(es - ea, 0.0)
end
