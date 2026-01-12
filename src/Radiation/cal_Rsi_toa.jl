"""
  cal_Rsi_toa(lat, J)

Estimate daily extraterrestrial radiation in MJ m-2 day-1.

# Arguments
- `lat`: Latitude in degrees.
- `J`: Day of the year.

# Returns
- Extraterrestrial radiation in `MJ m-2 d-1`.
"""
function cal_Rsi_toa(lat=0, J::Integer=1)
  dr = 1 + 0.033 * cos(π * J / 182.5) # Allen, Eq. 23
  σ = 0.409 * sin(π * J / 182.5 - 1.39) # Allen, Eq. 24

  ws = HourAngleSunSet(lat, J)
  # 0.0820 MJ m-2 min-1
  # 24 * 60 * 0.082 = 118.08
  lat = deg2rad(lat)
  Rsi_toa = 118.08 * dr / π * (ws * sin(lat) * sin(σ) + cos(lat) * cos(σ) * sin(ws)) # Allen, Eq. 21
  Rsi_toa
  # max(MJ2W(Rsi_toa), 0)
end


"""
    cal_Rsi_toa_hour(hour_beg=12, hour_end=hour_beg + 1; lat=0, J::Integer=1)

In MJ m-2 per unit time.
"""
function cal_Rsi_toa_hour(hour_beg=12, hour_end=hour_beg + 1; lat=0, J::Integer=1)
  dr = 1 + 0.033 * cos(π * J / 182.5) # Allen, Eq. 23
  σ = 0.409 * sin(π * J / 182.5 - 1.39) # Allen, Eq. 24

  lat = deg2rad(lat)
  coef = 24 * 60 * 0.082 # [1440 min] × [0.082 MJ min-1]

  w_beg = (hour_beg - 12) / 24 * 2π
  w_end = (hour_end - 12) / 24 * 2π

  ws = HourAngleSunSet(lat, J)
  w_beg = clamp(w_beg, -ws, ws)
  w_end = clamp(w_end, -ws, ws)

  # MJ m-2
  Rsi_toa = coef * dr / 2π * (
    (w_end - w_beg) * sin(lat) * sin(σ) + cos(lat) * cos(σ) * (sin(w_end) - sin(w_beg))) # Allen, Eq. 28
  max(Rsi_toa, 0)
  # # MJ m-2
  # max(Rsi_toa*1e6/3600, 0)
end


function cal_Rsi_toa_hour(dates::Vector{DateTime}; lat=0)
  map(date -> begin
      date_beg = date
      date_end = date + Hour(1)
      hour_beg = hour(date_beg)
      hour_end = hour(date_end)
      J = dayofyear(date)
      Rs = cal_Rsi_toa_hour(hour_beg, hour_end; lat, J)
      Rs * 1e6 / 3600 # hourly average
    end, dates)
end

"""
    cal_Rsi_toa_inst(dates::Vector{DateTime}; lat=0)

[W m-2]
"""
function cal_Rsi_toa_inst(dates::Vector{DateTime}; lat=0)
  Rg = 1367.0 # [W m-2]
  map(date -> begin
      J = dayofyear(date)
      r_v = 1 + 0.033 * cos(π * J / 182.5) # Allen, Eq. 23, r_v = (d_0/d)^2
      h = angle_SunElevation(lat, date)
      sin_h = max(sin(h), 0.0)
      Rg * r_v * sin_h
    end, dates)
end


export cal_Rsi_toa, cal_Rsi_toa_hour, cal_Rsi_toa_inst
