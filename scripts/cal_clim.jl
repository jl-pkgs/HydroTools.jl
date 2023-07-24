function cal_climIndex(f; outdir=".", overwrite=false)
  prefix = str_extract(basename(f), ".*(?=_\\d{4})")
  fout = "$outdir/$(prefix)_climIndex.nc"
  
  if isfile(fout) && !overwrite
    println("File exists: $fout")
    return
  end
  
  dates = nc_date(f)
  data = nc_read(f)

  years = year.(dates)
  grps = unique(years)
  grps = grps[grps.<=2014]

  ## 生成reference period
  period_ref = [1961, 1990]
  inds_ref = period_ref[1] .<= years .<= period_ref[2]
  data_ref = selectdim(data, 3, inds_ref)
  q95 = mapslices(x -> NanQuantile(x; probs=[0.95]), data_ref, dims=3)[:, :, 1] |> x -> Float32.(x)
  # q95 = Float32.(q95)

  ## 计算
  nlon, nlat = 140, 80
  varnames = ["cdd", "r20mm", "rx5day", "r95ptot", "prcptot"]
  nvar = length(varnames)
  res = zeros(Float32, nlon, nlat, length(grps), nvar)

  @time @inbounds for k in eachindex(grps)
    year = grps[k]
    mod(k, 10) == 0 && println("year = $year")
    ind = findall(years .== year)
    @views x = data[:, :, ind]
    x5 = zeros(Float32, size(x, 3) - 4)

    for i = 1:nlon, j = 1:nlat
      res[i, j, k, :] = index_P!(x5, x[i, j, :], q95[i, j])
    end
  end

  dims = [
    ncvar_dim(f)[1:2]...
    # NcDim("lon", xx, Dict("longname" => "Longitude", "units" => "degrees east"))
    # NcDim("lat", yy, Dict("longname" => "Latitude", "units" => "degrees north"))
    NcDim("year", grps)
    NcDim("index", varnames)
    # nc_dim(nc, "time")
  ]

  println("Writing ...")
  band = nc_bands(f)[1]
  @time nc_write(fout, band, res, dims)
end

# using BenchmarkTools
# x = rand(Float64, 365)
# x5 = zeros(Float64, 361)
# @benchmark index_P(x, 0.95)
# @benchmark index_P!(x5, x, 0.95)
# @profview index_P(x, 0.95)
