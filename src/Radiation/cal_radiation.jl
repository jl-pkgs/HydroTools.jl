export cal_Rn, cal_Rsi, cal_Rsi_toa, 
  cal_Rln, cal_Rln_out, cal_Rli, cal_Rln_yang2019
export blackbody


"""
blackbody radiation

# Arguments
- `Ta`: air temperature, in degC
- `ϵ`: emissivity, default is 1

# Return
- longwave radiation, in `W m-2`
"""
function blackbody(ϵ::T, Ta::T) where {T<:Real}
  ϵ * σ * (Ta + K0)^4
end


# 辐射用的单位都是W
"""
    cal_Rn(Rs, Rln, Tavg, albedo, emiss)

Calculate net radiation (Rn) using input parameters.

# Arguments

- `Rs`   : Shortwave inward solar radiation, W m-2.
- `Rln`  : Longwave inward solar radiation, W m-2.
- `Tavg` : Average temperature in Celsius.
- `α`    : Shortwave albedo.
- `ϵ`    : Longwave emissivity.

# Returns
- `Rn   `: Net radiation, W m-2.
"""
function cal_Rn(Rs::T, Rln::T, Tavg::T, α::Float64, ϵ::Float64) where {T<:Real}
  RLout = Stefan * (Tavg + 273.15)^4 / 0.0864 # convert from MJ m-2 d-1 to W/m2

  Rnl = ϵ * (Rln - RLout)
  Rns = (1.0 - α) * Rs

  Rn = Rns + Rnl
  Rn
  # Rn, Rns, Rnl
  # max(Rn, T(0.0))
end

cal_Rln_out(T::FT, ϵ=1.0) where {FT<:Real} = ϵ * σ * (T + K0)^4

"""
  cal_Rsi_toa(lat, J)

Estimate daily extraterrestrial radiation in MJ m-2 day-1.

# Arguments
- `lat`: Latitude in degrees.
- `J`: Day of the year.

# Returns
- Extraterrestrial radiation in `W m-2`.
"""
function cal_Rsi_toa(lat=0, J::Integer=1)
  dr = 1 + 0.033 * cos(π * J / 182.5) # Allen, Eq. 23
  σ = 0.409 * sin(π * J / 182.5 - 1.39) # Allen, Eq. 24

  ws = HourAngleSunSet(lat, J)
  # 24 * 60 * 0.082 = 118.08
  lat = deg2rad(lat)
  Rsi_toa = 118.08 * dr / π * (ws * sin(lat) * sin(σ) + cos(lat) * cos(σ) * sin(ws)) # Allen, Eq. 21
  max(MJ2W(Rsi_toa), 0)
end


"""
  cal_Rsi(lat, J, ssd=nothing, cld=nothing, Z=0, a=0.25, b=0.5)

Daily inward shortwave solar radiation at crop surface in W m-2 by
providing sunshine duration (SSD) in hours or cloud cover in fraction.

# Arguments
- `lat`: Latitude in degrees.
- `J`: Day of the year.
- `ssd`: Sunshine duration in hours. If `ssd` is `nothing`, `Rsi` is the
  clear-sky solar radiation. Default is `nothing`.
- `cld`: Cloud cover in fraction. If provided it would be directly used to
  calculate solar radiation rather than SSD and parameter a and b. Default is
  `nothing`.
- `Z`: Elevation in meters. Default is 0.
- `a`: Coefficient of the Angstrom formula. Determine the relationship between
  ssd and radiation. Default is 0.25.
- `b`: Coefficient of the Angstrom formula. Default is 0.50.

# Returns
- A tuple of solar radiation at crop surface in W m-2:
  - `Rsi`: Surface downward shortwave radiation.
  - `Rsi_o`: Clear-sky surface downward shortwave radiation.
  - `Rsi_toa`: Extraterrestrial radiation.
"""
function cal_Rsi(lat=0, J::Integer=1, ssd=nothing, cld=nothing; 
  Z=0, a=0.25, b=0.5)

  Rsi_toa = cal_Rsi_toa(lat, J)

  if cld !== nothing
    nN = (1 - cld)
  else
    N = SunshineDuration(lat, J)
    ssd = ssd === nothing ? N : ssd
    nN = ssd / N
  end

  Rsi = (a + b * nN) * Rsi_toa
  Rsi_o = (0.75 + 2 * Z / 1e5) * Rsi_toa

  Rsi, Rsi_o, Rsi_toa
end


"""
    cal_Rln(Tmax, Tmin, ea, cld)
    cal_Rln(Tmax, Tmin, ea, Rsi, Rsi_o)

Net outgoing longwave radiation.

# Arguments
- `Tmax`: Daily maximum air temperature at 2m height in degrees Celsius.
- `Tmin`: Daily minimum air temperature at 2m height in degrees Celsius.
- `ea`: Actual vapor pressure in kPa. Can be estimated by maximum or minimum air
  temperature and mean relative humidity.
- `Rsi`: Incoming shortwave radiation at crop surface in MJ m-2 day-1. Default
  is `nothing`.
- `Rsi_o`: Clear sky incoming shortwave radiation, i.e., extraterrestrial
  radiation multiply by clear sky transmissivity (i.e., a + b, a and b are
  coefficients of Angstrom formula. Normally 0.75) in MJ m-2 day-1. Default is
  `nothing`.
- `cld`: Cloud cover in fraction. If provided it would be directly used to
  calculate net outgoing longwave radiation than Rso. Default is `nothing`.

# Returns
- Net outgoing longwave radiation in MJ m-2 day-1.
"""
function cal_Rln(Tmax, Tmin, ea, cld)
  cld = isnan(cld) ? 1.0 : cld
  return (4.093e-9 * (((Tmax + 273.15)^4 + (Tmin + 273.15)^4) / 2)) *
         (0.34 - (0.14 * sqrt(ea))) *
         (1.35 * (1.0 - cld) - 0.35)
end

function cal_Rln(Tmax, Tmin, ea, Rsi, Rsi_o)
  cld = 1.0 - Rsi / Rsi_o
  cal_Rln(Tmax, Tmin, ea, cld)
end


"""
    cal_Rln_yang2019(Tₛ, Rsi, Rsi_toa; lat=30, ϵ::Float64=0.96, n1=2.52, n2=2.37, n3=0.035)

Calculate the longwave radiation based on the Yang et al. (2019) method.

# Arguments
- `Ts`: Land surface temperature.
- `Rsi`: Surface incoming short-wave radiation.
- `Rsi_toa`: Incoming short-wave radiation at the top of the atmosphere.
- `lat`: Latitude (in degrees). Default is 30.
- `ϵ`: Emissivity. Default is 0.96.
- `n1`, `n2`, `n3`: Parameters for `ΔT`. See Yang 2019 for details.

# References
1. Yang, Y., & Roderick, M. L. (2019). Radiation, surface temperature and
   evaporation over wet surfaces. Quarterly Journal of the Royal Meteorological
   Society, 145(720), 1118–1129. https://doi.org/10.1002/qj.3481
2. Tu, Z., & Yang, Y. (2022). On the Estimation of Potential Evaporation Under
   Wet and Dry Conditions. Water Resources Research, 58(4).
   https://doi.org/10.1029/2021wr031486
"""
function cal_Rln_yang2019(Ts, Rsi, Rsi_toa;
  lat=30, ϵ::Float64=0.96, n1=2.52, n2=2.37, n3=0.035)

  tau = Rsi / Rsi_toa
  ΔT = n1 * exp(n2 * tau) + n3 * abs(lat)
  return ϵ * σ * ((Ts - ΔT)^4 - Ts^4)
end


"""
  cal_Rli(Ta, ea=0, s=1, method="MAR")

Estimate Incoming longwave radiation. Not included in FAO56 paper but added for
convenience.

# Arguments
- `Ta`: Near surface air temperature in degrees Celsius.
- `ea`: Near surface actual vapour pressure in kPa. Default is 0.
- `s`: Cloud air emissivity, the ratio between actual incoming shortwave
  radiation and clear sky incoming shortwave radiation. Default is 1.
- `method`: The method to estimate the air emissivity. Must be one of 'MAR',
  'SWI', 'IJ', 'BRU', 'SAT', 'KON'. Default is 'MAR'.

# Returns
- Incoming longwave radiation in W/m².

# References
1. Satterlund, D. R. (1979), An improved equation for estimating long-wave
   radiation from the atmosphere, Water Resour. Res., 15( 6), 1649– 1650,
   doi:10.1029/WR015i006p01649.
2. Sedlar, J., & Hock, R. (2009). Testing longwave radiation parameterizations
   under clear and overcast skies at Storglaciären, Sweden. The Cryosphere,
   3(1), 75-84.
"""
function cal_Rli(Ta, ea=0; s=1, method="MAR")
  ea *= 10 # to hPa
  Ta += 273.15
  if method == "MAR"
    ep_ac = 0.5893 + 5.351e-2 * sqrt(ea / 10)
  elseif method == "SWI"
    ep_ac = 9.294e-6 * Ta^2
  elseif method == "IJ"
    ep_ac = 1 - 0.261 * exp(-7.77e-4 * (273 - Ta)^2)
  elseif method == "BRU"
    ep_ac = 1.24 * (ea / Ta)^(1 / 7)
  elseif method == "SAT"
    ep_ac = 1.08 * (1 - exp(-ea^(Ta / 2016)))
  elseif method == "KON"
    a = 0.4393
    b = 7
    ep_ac = 0.23 + a * (ea / 10 / Ta)^(1 / b)
  end
  ep_a = 1 - s + s * ep_ac
  return ep_a * σ * Ta^4
end
