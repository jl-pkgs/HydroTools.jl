#' @param z.wind Height where measure the wind speed `[m]`. Default 10m.
#' @param Uz wind speed at height `z.wind`
#' @param U2 wind speed at 2m
function cal_U2(Uz::T, z_wind=10.0) where {T<:Real}
  z_wind == 2.0 ? Uz : Uz * 4.87 / log(67.8 * z_wind - 5.42)
end

function cal_lambda(Tair::T) where {T<:Real}
  (2500.0 - Tair * 2.2) / 1000.0 * u"MJ / kg"
end

function cal_slope(Tair::T) where {T<:Real}
  4098.0 * (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) / (Tair + 237.3)^2 * u"kPa / K"
end


export cal_lambda, cal_slope
