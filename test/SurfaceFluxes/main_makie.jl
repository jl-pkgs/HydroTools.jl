using MakieLayers, GLMakie, Printf

function plot_soil(Tsoil)
  t = 0.5:0.5:24
  x = 1:nday
  y = t

  xticks = 0:30:nday |> format_ticks
  yticks = 0:3:24 |> format_ticks
  # ticks = 0:0.1:0.5 |> format_ticks

  layers = [1:6; [8, 9, 10]]
  fig = Figure(; size=(1400, 800))

  titles = [@sprintf("layer=%d, depth=%.3f", i, -soil.z[i]) for i in layers]
  imagesc!(fig, x, y, Tsoil[:, :, layers];
    colorbar=(; width=15),
    axis=(; xticks, yticks), titles, byrow=true)
  labs!(fig, xlabel="Day", ylabel="Hour"; fontsize=16, height=8)
  fig
end

# begin
#   fig = Figure(; size=(1400, 600))
#   ax = Axis(fig[1, 1]; xlabel="day", ylabel="hour")
#   lines!(ax, 1:10; label="1")
#   axislegend(ax)
#   fig
# end
