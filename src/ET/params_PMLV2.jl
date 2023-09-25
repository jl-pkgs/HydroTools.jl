using Parameters

# canopy height, in the order of `IGBP006` code
hc_raw = [10, 10, 10, 10, 10, 1, 1, 5, 5, 0.2, 1, 0.5, 10, 1, 0.01, 0.05, 0.01]

FT = Float64

parNames = ["Alpha", "Thelta", "m", "Am_25",
  "VPDmin", "VPDmax", "D0", "kQ", "kA",
  "S_sls", "fER0", 
  "hc"]

params = [
  0.01 0.10 0.06  # `Alpha` : initial slope of the light response curve to assimilation 
                  #           rate (i.e., quantum efficiency; μmol CO2 [μmol PAR]−1)
  0.01 0.07 0.04  # `Thelta`: initial slope of the CO2 response curve to assimilation 
                  #           rate (i.e., carboxyla-tion efficiency; μmol m−2 s−1 [μmol m−2 s−1]−1),
  2.00 100 10     # `m`     : stomatal conductance coefficient,
  5.00 120 50     # `Am_25` : carbon saturated rate of photosynthesis at 25 °C, μmol m−2 s−1
  0.65 1.5 0.9    # `VPDmin`: kPa
  3.50 6.5 4.0    # `VPDmin`: kPa
  0.50 2.0 0.7    # `D0`
  0.10 1.0 0.45   # `kQ`    : extinction coefficients for visible radiation
  0.50 0.9 0.7    # `kA`    : extinction coefficients for available energy
  0.01 1.0 0.1    # `S_sls` : 
  0.01 0.5 0.1    # `fER0`  : 
  0.01 20  1      # hc
  # 1 6 4           # `LAIref`
  # 6 14 10         # `frame`: 8-day moving window
]

theta0 = params[:, 3]
parRanges = params[:, 1:2]
param0 = list(parNames, theta0)

theta2par(theta) = list(parNames, theta)

export parRanges, parNames, theta0, param0, hc_raw
export param_PML, theta2par;
