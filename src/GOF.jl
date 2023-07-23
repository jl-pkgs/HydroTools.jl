import Statistics: mean, std, cor


valid_index(x::AbstractVector) = @.(!(ismissing(x) || isnan(x)))

function valid_index(obs::AbstractVector, sim::AbstractVector)
  inds = valid_index(obs) .&& valid_index(sim)

  obs = @view obs[inds]
  sim = @view sim[inds]
  obs, sim
end


"""
`of_KGE` Calculates Kling-Gupta Efficiency of simulated streamflow (Gupta et al,
2009). Ignores time steps with negative flow values.

# References

1. Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009).
Decomposition of the mean squared error and NSE performance criteria:
Implications for improving hydrological modelling. Journal of Hydrology,
377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
"""
function of_KGE(obs, sim, w=[1, 1, 1])
  obs, sim = valid_index(obs, sim)
  length(sim) <= 2 && (return -999.0)
  ## Check inputs and select timesteps
  
  ## calculate components
  c1 = cor(obs, sim)                    # r: linear correlation
  c2 = std(sim) / std(obs)               # alpha: ratio of standard deviations
  c3 = mean(sim) / mean(obs)             # beta: bias 

  ## calculate value
  1 - sqrt((w[1] * (c1 - 1))^2 + (w[2] * (c2 - 1))^2 + (w[3] * (c3 - 1))^2)    # weighted KGE
end


function of_NSE(obs, sim)
  obs, sim = valid_index(obs, sim)
  length(sim) <= 2 && (return -999.0)
  
  # %% Calculate metric
  top = sum((sim .- obs) .^ 2)
  bot = sum((obs .- mean(obs)) .^ 2)
  1 - (top / bot)
end


# both of low and high
function of_KGE_multi(obs, sim; min=0.01)
  sim[sim .< min] .= min
  of_KGE(obs, sim) + of_KGE(1.0 ./ obs, 1.0 ./ sim)
end

function of_NSE_multi(obs, sim; min=0.01)
  # obs[obs.<0.01] = 0.01
  sim[sim .< min] .= min
  of_NSE(obs, sim) + of_NSE(1.0 ./ obs, 1.0 ./ sim)
end


function GOF(obs::AbstractVector{T}, sim::AbstractVector{T}) where {T<:Real}
  obs, sim = valid_index(obs, sim)

  n = length(obs)
  n_valid = n
  e = sim - obs
  μ = mean(obs)

  bias = mean(e)
  bias_perc = bias / μ * 100
  RMSE = sqrt(sum(e .^ 2) / n)
  MAE = mean(abs.(e))


  if length(sim) <= 2 
    return (; KGE=-999.0, NSE=-999.0, R2=NaN, R=NaN, RMSE, MAE, bias, bias_perc, n_valid)
  end

  KGE = of_KGE(obs, sim)
  NSE = of_NSE(obs, sim)
  R = cor(obs, sim)
  R2 = R^2
  
  (; KGE, NSE, R2, R, RMSE, MAE, bias, bias_perc, n_valid)
end


export of_NSE, of_KGE, of_NSE_multi, of_KGE_multi, 
  GOF, valid_index
