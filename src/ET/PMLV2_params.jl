using Parameters
using FieldMetadata
import FieldMetadata: @bounds, bounds
export bounds, Param_PMLV2


abstract type AbstractETParam{FT<:AbstractFloat} end


"""
    struct Param_PMLV2{FT<:AbstractFloat} <: AbstractETParam{FT}

# Fields
$(TYPEDFIELDS)
"""
@bounds @with_kw mutable struct Param_PMLV2{FT<:AbstractFloat} <: AbstractETParam{FT}
  "initial slope of the light response curve to assimilation rate, (i.e., quantum efficiency; `μmol CO2 [μmol PAR]⁻¹`)`"
  Alpha ::FT = 0.06  | (0.01, 0.10)

  "initial slope of the CO2 response curve to assimilation rate, (i.e., carboxylation efficiency; `μmol m⁻² s⁻¹ [μmol m⁻² s⁻¹]⁻¹`)"
  Thelta::FT = 0.04  | (0.01, 0.07)

  "stomatal conductance coefficient"
  m     ::FT = 10.00 | (2.00, 100.00)

  "carbon saturated rate of photosynthesis at 25 °C, `μmol m⁻² s⁻¹`"
  Am_25 ::FT = 50.00 | (5.00, 120.00)
  
  "parameter to constrain `gc`, kPa"
  VPDmin::FT = 0.9   | (0.65, 1.5)
  "parameter to constrain `gc`, kPa"
  VPDmax::FT = 4.0   | (3.50, 6.5)

  ""
  D0    ::FT = 0.7   | (0.50, 2.0)
  "extinction coefficients for visible radiation"
  kQ    ::FT = 0.45  | (0.10, 1.0)
  "extinction coefficients for available energy"
  kA    ::FT = 0.70  | (0.50, 0.9)

  "Specific leaf storage, van Dijk, A.I.J.M, 2001, Eq2"
  S_sls ::FT = 0.1   | (0.01, 1.0)
  "Canopy cover fraction related parameter"
  fER0  ::FT = 0.1   | (0.01, 0.5)
  "canopy height, `[m]`"
  hc    ::FT = 1.0   | (0.01, 20.0)
  # LAIref::FT = 4.0   | (1.0, 6.0)      # 
  # frame::Integer = 10.0  | (6.0, 14.0) # 8-day moving window
end


function Base.collect(par::AbstractETParam)
  [getfield(par, f) for f in fieldnames(typeof(par))]
end

function get_bounds(par::AbstractETParam)
  hcat(map(collect, bounds(par))...) |> transpose |> collect
end

# canopy height, in the order of `IGBP006` code
hc_raw = [10, 10, 10, 10, 10, 1, 1, 5, 5, 0.2, 1, 0.5, 10, 1, 0.01, 0.05, 0.01]

FT = Float64

par0 = Param_PMLV2()
theta0 = collect(par0)
parNames = fieldnames(Param_PMLV2)
parRanges = get_bounds(par0)

param0 = list(parNames, theta0)

theta2par(theta) = list(parNames, theta)

export parRanges, parNames, theta0, param0, hc_raw
export param_PML, theta2par;
