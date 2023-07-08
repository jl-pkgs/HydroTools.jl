function cal_es(Tair::T)::T where {T<:Real}
  # FAO98 equation
  0.6108 * exp((17.27 * Tair) / (Tair + 237.3))
  # 0.61094 * exp((17.625 * Tair) / (Tair + 243.04))
end

function Tdew2RH(Tdew::T, Tair::T)::T where {T<:Real}
  ea::T = cal_es(Tdew)
  es::T = cal_es(Tair)
  ea / es * 100
end

function Tdew2VPD(Tdew::T, Tair::T)::T where {T<:Real}
  ea::T = cal_es(Tdew)
  es::T = cal_es(Tair)
  es - ea
end

function q2ea(q, Pa=atm)
  q * Pa/(epsilon + (1 - epsilon) * q)
end

function q2RH(q, Tair, Pa=atm)
  ea = q2ea(q, Pa)
  es = cal_es(Tair)
  ea / es * 100
end
