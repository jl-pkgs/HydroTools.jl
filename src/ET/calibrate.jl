function m_goal(df, theta; IGBPcode=nothing, of_gof=:NSE, verbose=false)
  # IGBPcode !== nothing && (par.hc = hc_raw[IGBPcode])
  IGBPcode !== nothing && (theta[end] = hc_raw[IGBPcode]) # the last one is hc
  par = list(parNames, theta)
  # @show theta
  
  dobs = df[!, [:GPP_obs, :ET_obs]]
  dsim = PMLV2_sites(df; par)

  ## 8-day, yearly, yearly_anom
  # 1. 8-day
  info_GPP = GOF(dobs.GPP_obs, dsim.GPP)
  info_ET = GOF(dobs.ET_obs, dsim.ET)

  if verbose
    indexes = [:KGE, :NSE, :R2, :RMSE, :bias, :bias_perc]
    info_GPP = info_GPP[indexes] |> round2
    info_ET  = info_ET[indexes] |> round2

    @show info_GPP
    @show info_ET
  end

  goal = (info_ET[of_gof] + info_GPP[of_gof]) / 2
  goal
end


## 最后一步，参数率定模块
function m_calib(df::AbstractDataFrame; IGBPcode=nothing, maxn=2500, kw...)
  theta, goal, flag = sceua(theta -> -m_goal(df, theta; IGBPcode, kw...),
    theta0, parRanges[:, 1], parRanges[:, 2]; maxn, kw...)
  theta, goal, flag
end


# rounded_data = NamedTuple((field => round(value) for (field, value) in data))

export m_goal, m_calib
