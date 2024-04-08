function get_md(J)
  date_origin = DateTime(2010, 1, 1, 0, 0, 0) + Day(J - 1)
  Dates.format(date_origin, "mm-dd")
end

function ssh2time(ssh)
  w_hour = ssh / 2
  w_sec = floor(Int, w_hour * 3600)
  time_begin = DateTime(2000, 1, 1, 12, 0, 0) - Dates.Second(w_sec)
  time_end = DateTime(2000, 1, 1, 12, 0, 0) + Dates.Second(w_sec)
  # time_begin, time_end
  time_begin, time_end = Dates.format(time_begin, "HH:MM:SS"),
  Dates.format(time_end, "HH:MM:SS")
  time_begin, time_end
end


"""
    SolarDeclinationAngle(J; to_deg=false)

# Arguments
- σ: Solar Declination Angle, 黄赤交角（太阳赤纬角）
- ϕ: `纬度`
- ω: `时角`
"""
function SolarDeclinationAngle(J::Integer; deg=false)
  σ = 0.409 * sin(2 * pi / 365 * J .- 1.39) # Allen, Eq. 24
  deg ? rad2deg(σ) : σ
end

"""
  HourAngleSunSet(lat, J)

Calculate the sunset hour angle.

# Arguments
- `lat`: Latitude in degrees.
- `J`: Day of the year.

# Returns
- Sunset hour angle in radian
"""
function HourAngleSunSet(lat=0, J::Integer=1; deg=false)
  lat = deg2rad(lat)
  σ = SolarDeclinationAngle(J)

  tmp = clamp(-tan(lat) * tan(σ), -1, 1)
  ws = acos(tmp) # Eq. 25
  deg ? rad2deg(ws) : ws
end

function SunshineDuration(lat=0, J::Integer=1)
  w = HourAngleSunSet(lat, J)
  w_hour = w / pi * 12 # 距离中午的时间
  w_hour * 2 # ssh
end

export HourAngleSunSet, SunshineDuration, ssh2time
