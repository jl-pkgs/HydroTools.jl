deg2rad(x) = x / 180 * pi
rad2deg(x) = x / pi * 180


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
    get_σ(J; to_deg=false)


- `黄赤交角` : σ, Solar Declination Angle

- `纬度`: ϕ

- `时角`: ω
"""
function get_σ(J; to_deg=false)
  σ = 0.409 .* sin.(2 * pi / 365 * J .- 1.39) # in [rad]
  if to_deg
    σ = rad2deg.(σ)
  end
  σ
end

function get_ssh(ϕ, J)
  σ = get_σ.(J)
  # rad2deg(σ)
  tmp = -tan.(ϕ) .* tan.(σ)
  clamp!(tmp, -1, 1)
  # constrain in the range of [-1, 1]
  w = acos.(tmp)
  w_hour = w / pi * 12 # 距离中午的时间
  w_hour * 2 # ssh
end
