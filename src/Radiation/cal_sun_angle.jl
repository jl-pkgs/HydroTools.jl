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
    angle_SolarDeclination(J; to_deg=false)

# Arguments
- σ: Solar Declination Angle, 黄赤交角（太阳赤纬角）
- ϕ: `纬度`
- ω: `时角`
"""
function angle_SolarDeclination(J::Integer; deg=false)
  σ = 0.409 * sin(2pi / 365 * J .- 1.39) # Allen, Eq. 24
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
@fastmath function HourAngleSunSet(lat=0, J::Integer=1; deg=false)
  lat = deg2rad(lat)
  σ = angle_SolarDeclination(J)

  tmp = clamp(-tan(lat) * tan(σ), -1, 1)
  ws = acos(tmp) # Eq. 25
  deg ? rad2deg(ws) : ws
end

function SunshineDuration(lat=0, J::Integer=1)
  w = HourAngleSunSet(lat, J)
  w_hour = w / pi * 12 # 距离中午的时间
  w_hour * 2 # ssh
end


"""
任一点，任一时刻的太阳高度角
"""
function angle_SunElevation(lat::Real, time_local::DateTime)
  ψ = deg2rad(lat)  # 纬度转弧度
  J = dayofyear(time_local)
  σ = angle_SolarDeclination(J; deg=false)
  dh = hour(time_local) + minute(time_local) / 60 - 12.0
  ω = deg2rad(dh * 15) # 时角，上午为负，下午为正
  # @show ψ, J, σ, ω
  sinH = cos(ψ) * cos(σ) * cos(ω) + sin(ψ) * sin(σ)
  return asin(sinH)
end


export HourAngleSunSet, SunshineDuration, ssh2time,
  angle_SunElevation
