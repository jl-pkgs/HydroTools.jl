# calibrate for PMLV2 parameter
using DataFrames

function to_df(res::output_PML{T}) where {T<:Real}
  data, names = to_mat(res)
  DataFrame(data, names)
end

# 一个站点的计算
function PMLV2(Prcp::AbstractVector{T}, Tavg::AbstractVector{T},
  Rs::AbstractVector{T}, Rn::AbstractVector{T},
  VPD::AbstractVector{T}, U2::AbstractVector{T}, LAI::AbstractVector{T},
  Pa::AbstractVector{T}, Ca::AbstractVector{T};
  par=param0, frame=3,
  res::Union{Nothing,output_PML}=nothing) where {T<:Real}

  n = length(Prcp)
  fields = fieldnames(interm_PML)[2:end-2] # skip ET, fval_soil and Es
  r = interm_PML{T}()
  res === nothing && (res = output_PML{T}(; n))

  for t = 1:n
    PMLV2(Prcp[t], Tavg[t], Rs[t], Rn[t], VPD[t], U2[t], LAI[t], Pa[t], Ca[t]; par, r)
    res[t, fields] = r
  end

  res.fval_soil = movmean2(Prcp, frame, 0) ./ movmean2(res.Es_eq, frame, 0)
  clamp!(res.fval_soil, 0, 1)

  res.Es .= res.fval_soil .* res.Es_eq
  res.ET .= res.Ec .+ res.Ei .+ res.Es
  res
end



function PMLV2(d::AbstractDataFrame; par=param0, kw...)
  PMLV2(d.Prcp, d.Tavg, d.Rs, d.Rn,
    d.VPD, d.U2, d.LAI,
    d.Pa, d.Ca; par, kw...) |> to_df
end


# 相同植被类型多个站点一起的计算
function PMLV2_sites(df::AbstractDataFrame; par=param0, kw...)
  sites = df.site
  grps = unique(sites)

  # 这里有优化空间，但程序会写的很复杂
  res = []
  for grp in grps
    inds = sites .== grp
    d = df[inds, :]
    r = PMLV2(d; par, kw...)
    push!(res, r)
  end  
  vcat(res...) 
end


function m_goal(df, theta; IGBPcode=nothing, of_gof=:NSE, verbose=false)
  par = list(parNames, theta)
  IGBPcode !== nothing && (par.hc = hc_raw[IGBPcode])

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
  theta, goal, flag = sceua(theta -> -m_goal(df, theta; kw...),
    theta0, parRanges[:, 1], parRanges[:, 2]; maxn, kw...)
  theta, goal, flag
end


round2(x::NamedTuple, digits=3; kw...) = map(val -> round(val; digits=digits), x)

export round2;
# rounded_data = NamedTuple((field => round(value) for (field, value) in data))


export PMLV2_sites, 
  m_goal, m_calib
