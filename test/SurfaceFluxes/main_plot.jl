using RTableTools
using Printf
using Plots
import Plots: mm
using HydroTools, Test
using HydroTools.SurfaceFluxes
# pyplot()

gr(framestyle=:box)
d_mat = fread("$(@__DIR__)/dat_Fluxes_day1_Ts=25.csv")

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
# function plot_soil(day)
#   z = Tsoil[day, :, :] .- K0
#   heatmap(t, -soil.z, z', title="day=$day",
#     yflip=true, yaxis=:log10, c=:jet)
# end

function plot_soil_layer(layer)
  z = Tsoil[:, :, layer] .- K0
  x = t
  y = 1:nday
  title = @sprintf("layer=%d, depth=%.3f", layer, -soil.z[layer])
  heatmap(y, x, z'; title, yflip=true, c=:jet)
end


Figure3_SOIL = plot(
  [plot_soil_layer(i) for i = [1:6; [8, 9, 10]]]...,
  size=(900, 700),
  ylabel="hour", xlabel="days",
  left_margin=0mm, right_margin=0mm, top_margin=1mm, bottom_margin=1mm,
  subplot_spacing=0mm
)
# fig = plot_soil(Tsoil)
# save("sp07_soil_temperature_Ts0=$Ts0.png", fig)

# check with MATLAB
Figure1 = plot(
  plot_comp_mat("Rn"),
  plot_comp_mat("H"),
  plot_comp_mat("LE"),
  plot_comp_mat("G"),
  plot_comp_mat("g_ac"),
  size=(1000, 600)
)


Figure2 = plot(
  plot_vary_day("Rn"),
  plot_vary_day("H"),
  plot_vary_day("LE"),
  plot_vary_day("G"),
  plot_vary_day("g_ac"),
  size=(1000, 600)
)
