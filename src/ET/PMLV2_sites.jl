# 相同植被类型多个站点一起的计算
function PMLV2_sites(df::AbstractDataFrame; par::Param_PMLV2=Param_PMLV2(), kw...)
  sites = df.site
  grps = unique(sites)

  # ! 这里有优化空间，但程序会写的很复杂
  # 每个站点单独率定
  res = []
  for grp in grps
    inds = sites .== grp
    d = df[inds, :]
    r = PMLV2(d; par, kw...)
    push!(res, r)
  end
  vcat(res...)
end


export PMLV2_sites
