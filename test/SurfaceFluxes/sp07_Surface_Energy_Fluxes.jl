using HydroTools, Test
using HydroTools.SurfaceFluxes
using RTableTools
# includet("main_plot.jl")
includet("main_makie.jl")

function model(soil, flux; nday=30)
  dt = 1800 # second
  ntime = Int(86400 / dt) # 48

  Rn = zeros(nday, ntime)
  H = zeros(nday, ntime)
  LE = zeros(nday, ntime)
  G = zeros(nday, ntime)
  g_ac = zeros(nday, ntime)
  Tsoil = zeros(nday, ntime, length(soil.dz))

  for i = 1:nday
    for j = 1:48
      hour = j * dt / 86400 * 24
      # printstyled("[j=$j] hour=$hour\n", color=:green, bold=true)
      # MET
      Ta = gen_Ta(hour)  # Air temperature (C)
      es, d_es = satvap(Ta) # Pa
      ea = es * RH / 100
      met = Met(Ta, ea, Pa, z; rain=0, snow=0, u=3.0)

      ## Radiation and Canopy
      Θ = Ta + K0
      Rln_in = (0.398e-05 * Θ^2.148) * σ * Θ^4
      rad = Radiation(day, hour, lat; Rln_in)

      coszen = cal_coszen(day, hour, lat)
      can = Canopy{Float64}(; LAI=5.0, coszen)

      surface_fluxes!(flux, met, rad, can, soil; param)
      Rn[i, j] = flux.Rn
      H[i, j] = flux.H
      LE[i, j] = flux.LE
      G[i, j] = flux.G_soil + flux.G_snow
      g_ac[i, j] = flux.g_ac
      Tsoil[i, j, :] .= soil.Tsoil
    end
  end
  return (; Rn, H, LE, G, g_ac=g_ac * 100, Tsoil)
end

begin
  z = 30.0 # Reference height (m)
  hc = 20.0 # canopy height
  d = 0.67 * hc # Zero-plane displacement height (m)
  z0m = 0.13 * hc # Roughness length for momentum (m)
  z0c = 0.10 * z0m # Roughness length for scalars (m)
  param = (; z, d, z0m, z0c)

  ## 固定值
  day = 182
  lat = 40.0 * pi / 180  # Latitude (degrees -> radians) for solar radiation

  z = 30.0
  Tmean = 25.0
  Pa = atm * 1e3
  RH = 70.0  # Relative humidity (%)
  
  ## Flux, 连续型的变量，全部不需管
  Ts = Tmean + K0
  es, d_es = satvap(Ts - K0)
  flux = Flux{Float64}(; θ_surf=Ts, e_surf=es)

  ## Soil
  dz = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058]
  soil = Soil(dz)
  Ts0 = 5.0
  init_soil!(soil; Ts=Ts0) # 每次只更新`Tsoil`

  t = 0.5:0.5:24
  nday = 365
  days = 1:nday
  
  R = model(soil, flux; nday)
  Tsoil = R[:Tsoil]
  # fig = plot_soil(Tsoil)
  # save("sp07_soil_temperature_Ts0=$Ts0.png", fig)
end

d_mat = fread("$(@__DIR__)/dat_Fluxes_day1_Ts=25.csv")

maximum(abs.(d_mat[:, :Rn] - R[:Rn][:])) <= 1e-3
maximum(abs.(d_mat[:, :H] - R[:H][:])) <= 0.002
maximum(abs.(d_mat[:, :LE] - R[:LE][:])) <= 1e-3
maximum(abs.(d_mat[:, :G] - R[:G][:]))  <= 1e-3
maximum(abs.(d_mat[:, :g_ac] - R[:g_ac][:])) <= 0.006

using MAT
f = "C:/Users/hydro/github/CUG-hydro/bonanmodeling/sp_07_Surface Energy Fluxes/Tsoil_5deg.mat"
Tsoil2 = matread(f)["Tsoil"]

## 计算存在误差，但没有本质差别
fig = plot_soil(Tsoil .- K0)
fig = plot_soil(Tsoil2 .- K0)
# fig = plot_soil(Tsoil2 .- Tsoil)


# # check with MATLAB
# Figure1 = plot(
#   plot_comp_mat("Rn"),
#   plot_comp_mat("H"),
#   plot_comp_mat("LE"),
#   plot_comp_mat("G"),
#   plot_comp_mat("g_ac"),
#   size=(1000, 600)
# )

# Figure2 = plot(
#   plot_vary_day("Rn"),
#   plot_vary_day("H"),
#   plot_vary_day("LE"),
#   plot_vary_day("G"),
#   plot_vary_day("g_ac"),
#   size=(1000, 600)
# )
