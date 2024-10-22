using Plots
using RTableTools
using Printf
gr(framestyle=:box)

d_mat = fread("$(@__DIR__)/dat_Surface_Energy_Fluxes.csv")

function plot_comp_mat(var; day=1, inds=1:48)
  t = 0.5:0.5:24
  _ticks = 0:3:24
  xticks = _ticks
  p = plot(; xticks, title=string(var))
  _var = Symbol(var)

  x_jl = R[_var][day, :]
  plot!(p, t[inds], x_jl[inds], label="Julia")
  plot!(p, t[inds], d_mat[inds, _var], label="MATLAB")
  p
end

function plot_vary_day(var; inds=1:48, days=[1, 2, 3, 4, 5, 10, 30])
  x = R[Symbol(var)]
  t = 0.5:0.5:24
  _ticks = 0:3:24
  xticks = _ticks
  p = plot(; xticks, title=string(var))

  for day = days
    _x = x[day, :]
    plot!(p, t[inds], _x[inds], label="day = $day")
  end
  p
end

# Tsoil
function plot_soil(day)
  z = Tsoil[day, :, :] .- K0
  heatmap(t, -soil.z, z', title="day=$day",
    yflip=true, yaxis=:log10, c=:jet)
end

function plot_soil_layer(layer)
  z = Tsoil[:, :, layer] .- K0
  x = t
  y = 1:nday
  title = @sprintf("layer=%d, depth=%.3f", layer, -soil.z[layer])
  heatmap(y, x, z'; title, yflip=true, c=:jet)
end
