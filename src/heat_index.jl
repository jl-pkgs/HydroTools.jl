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
"""
function heat_index(Tair::T, RH::T) where {T<:Real}
  Tair = C2F(Tair) # Fdeg

  if (Tair <= 40)
    HI = Tair
  else
    alpha = 61 + ((Tair - 68) * 1.2) + (RH * 0.094)
    HI = 0.5 * (alpha + Tair)
    if (HI > 79)
      HI = -42.379 + 2.04901523 * Tair + 10.14333127 * RH -
           0.22475541 * Tair * RH - 6.83783 * 10^-3 * Tair^2 -
           5.481717 * 10^-2 * RH^2 + 1.22874 * 10^-3 * Tair^2 * RH +
           8.5282 * 10^-4 * Tair * RH^2 - 1.99 * 10^-6 * Tair^2 * RH^2
      if (RH <= 13 && Tair >= 80 && Tair <= 112)
        adj1 = (13 - RH) / 4.0
        adj2 = sqrt((17 - abs(Tair - 95)) / 17.0)
        tol_adj = adj1 * adj2
        HI = HI - tol_adj
      elseif (RH > 85 && Tair >= 80 && Tair <= 87)
        adj1 = (RH - 85.0) / 10.0
        adj2 = (87.0 - Tair) / 5.0
        tol_adj = adj1 * adj2
        HI = HI + tol_adj
      end
    end
  end
  F2C(HI) # Cdeg
end
