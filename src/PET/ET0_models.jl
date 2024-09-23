# include("constant.jl")
# include("humidity.jl")


"""
    ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...)

Equilibrium evaporation, `ET0_eq = Δ / (Δ + γ) * Rn`

# Arguments
- `Rn`   : net radiation (W m-2)
- `Tair` : 2m air temperature (degC)
- `VPD`  : vapor pressure deficit (kPa)
- `Pa`   : surface air pressure (kPa)

# Returns
- `λ` : latent heat of vaporization (MJ kg-1)
- `Δ`  : slope of the saturation vapor pressure curve (kPa degC-1)
- `γ`  : psychrometric constant (kPa degC-1)
- `Eeq`    : equilibrium evaporation rate (mm day-1)

# Optional keyword arguments
- `args...` : additional arguments (not used in this function)

# Examples
```julia
julia> ET0_eq(200.0, 20.0)
(2.4536, 0.14474018811241365, 0.06725596686893338, 4.808405926729265)
```
"""
function ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...) where {T<:Real}
  # T = eltype(Rn)
  λ::T = cal_lambda(Tair) # MJ kg-1
  Δ::T = cal_slope(Tair) # kPa degC-1
  γ::T = Cp * Pa / (ϵ * λ) # kPa degC-1
  Eeq::T = Δ / (Δ + γ) * Rn |> x -> W2mm(x; λ)
  λ, Δ, γ, Eeq
end

"""
    ET0_Penman48(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2)

# Examples
```julia
# ET_water = ET0_Penman48(Rn, Tavg, VPD, U2, Pa)
ET0_Penman48(200., 20., 2., 2.)
```
"""
function ET0_Penman48(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2) where {T<:Real}
  # T = eltype(Rn)
  λ, Δ, γ, Eeq = ET0_eq(Rn, Tair, Pa)
  U2::T = cal_U2(wind, z_wind)
  # rho_a * Cp * dT / rH (MJ m-2 s-1)
  # rho_a ≈ 1.225 kg/m3
  # rho_a * Cp / rH = f(U2)
  # `f(U2) = 2.6 * (1 + 0.54U2)` is equivalent to Shuttleworth1993
  Evp::T = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ
  PET::T = Eeq + Evp
  PET
end


"""
    ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false)

# Examples
```julia
ET0_Penman48(200., 20., 2., 2.)
ET0_FAO98(200.0, 20.0, 2.0, 2.0)
```
"""
function ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false) where {T<:Real}
  # T = eltype(Rn)
  U2 = cal_U2(wind, z_wind)

  if tall_crop
    p1 = T(1600.0)
    p2 = T(0.38)
  else
    p1 = T(900.0)
    p2 = T(0.34)
  end
  λ, Δ, γ, Eeq = ET0_eq(Rn, Tair, Pa)

  Eeq::T = Δ / (Δ + (γ * (1.0 + p2 * U2))) * Rn |> x -> W2mm(x; λ)
  Evp::T = γ * p1 / (Tair + 273.15) * U2 * VPD / (Δ + (γ * (1.0 + p2 * U2)))

  PET::T = Eeq + Evp
  PET
end
