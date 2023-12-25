"""
    heat_index(Tair::T, RH::T) where {T<:Real}

Calculate the heat index given the air temperature `Tair` (in degrees Celsius) and the relative humidity `RH` (in percent).

# Arguments
- `Tair::T`: Air temperature in degrees Celsius.
- `RH::T`: Relative humidity in percent.

# Returns
- `HI::Float64`: Heat index in degrees Celsius.

# Examples
```julia-repl
julia> heat_index(30, 50)
31.049081444444305
```

# References

1. Steadman, R. G. “The Assessment of Sultriness. Part I: A Temperature-Humidity Index Based on Human Physiology and Clothing Science.” Journal of Applied Meteorology 18, no. 7 (July 1979): 861–73. https://doi.org/10.1175/1520-0450(1979)018<0861:TAOSPI>2.0.CO;2.

2. Steadman, Robert G. “A Universal Scale of Apparent Temperature.” Journal of Climate and Applied Meteorology 23, no. 12 (December 1984): 1674–87. https://doi.org/10.1175/1520-0450(1984)023<1674:AUSOAT>2.0.CO;2.

3. Rothfusz, Lans P. “The Heat Index ‘Equation’ (or, More Than You Ever Wanted to Know About Heat Index),” 1990.

4. Anderson, G. Brooke, Michelle L. Bell, and Roger D. Peng. “Methods to Calculate the Heat Index as an Exposure Metric in Environmental Health Research.” Environmental Health Perspectives 121, no. 10 (October 2013): 1111–19. https://doi.org/10.1289/ehp.1206273.
"""
function heat_index(Tair::T, RH::T) where {T<:Real}
  Tair = C2F(Tair) # Fdeg
  HI::T = Tair

  if (Tair <= 40) # about 4.44℃
    HI = Tair
  else
    alpha = 61 + ((Tair - 68) * 1.2) + (RH * 0.094)
    HI = 0.5 * (alpha + Tair)
    # HI = -10.3 + 1.1 * Tair + 0.047 * RH; Anderson2013, Figure3
    if (HI > 79)
      HI = _heat_index_Rothfusz1990(Tair, RH)
    end
  end
  F2C(HI) # Cdeg
end



function _heat_index_Rothfusz1990(F::T, RH::T) where {T<:Real}
  Tair = F
  # Note, Tair: Fdeg at here
  HI::T = -42.379 + 2.04901523 * Tair + 10.14333127 * RH -
       0.22475541 * Tair * RH - 6.83783 * 10^-3 * Tair^2 -
       5.481717 * 10^-2 * RH^2 + 1.22874 * 10^-3 * Tair^2 * RH +
       8.5282 * 10^-4 * Tair * RH^2 - 1.99 * 10^-6 * Tair^2 * RH^2
  if (RH <= 13 && 80 <= Tair <= 112)
    adj1 = (13 - RH) / 4.0
    adj2 = sqrt((17 - abs(Tair - 95)) / 17.0)
    tol_adj = adj1 * adj2
    HI = HI - tol_adj
  elseif (RH > 85 && 80 <= Tair <= 87)
    adj1 = (RH - 85.0) / 10.0
    adj2 = (87.0 - Tair) / 5.0
    tol_adj = adj1 * adj2
    HI = HI + tol_adj
  end
  HI
end


function wind_chill(Tair::T, U10::T) where {T<:Real}
  13.12 + 0.6215 * Tair - 11.37 * U10^0.16 + 0.3965 * Tair * U10^0.16
end

function universal_apparent_temperature(Tair::T, VPD::T, U10::T) where {T<:Real}
  -2.7 + 1.04 * Tair + 2 * VPD - 0.65 * U10 # Jianfeng, 2018, ncc, Eq. 3
end

function apparent_temperature(Tair::T, RH::T, U10::T) where {T<:Real}
  HI = heat_index(Tair, RH)
  WC = wind_chill(Tair, U10) 
  VPD = cal_es(Tair) * (1 - RH / 100)
  # 4.8 km h-1 = 1.333 m s-1
  
  if Tair <= 10 && U10 >= 1.333
    AP = WC
  elseif Tair >= 26
    AP = HI
  else
    AP = universal_apparent_temperature(Tair, VPD, U10)
  end
  AP
end

export _heat_index_Rothfusz1990, heat_index, wind_chill
