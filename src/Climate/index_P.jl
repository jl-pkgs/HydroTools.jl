function rollmax5!(x5, x::AbstractVector{T}) where {T<:Real}
  win = 5
  n = length(x)
  # x5 = zeros(T, n-win+1)
  @inbounds for i = 1:n-win+1
    x5[i] = max(x[i], x[i+1], x[i+2], x[i+3], x[i+4])
  end
  x5
end


function CDD(prcp::AbstractVector{T}) where {T<:Real}
  inds = findall(prcp .< T(1.0))
  flag = GroupFlag(inds)

  n = maximum(flag)
  len = zeros(Int, n)
  for i = 1:n
    len[i] = count(flag .== i)
  end
  # _, len = rle(grp)
  maximum(len)
end


function CWD(prcp::AbstractVector{T}) where {T<:Real}
  inds = findall(prcp .>= T(1.0))
  flag = GroupFlag(inds)
  
  n = maximum(flag)
  len = zeros(Int, n)
  for i = 1:n
    len[i] = count(flag .== i)
  end
  # _, len = rle(grp)
  maximum(len)
end

# function index_P(x, inds_ref::BitVector)
#   q95 = quantile(x[inds_ref], 0.95)
#   index_P(x, q95)
# end

function index_P!(x5::AbstractVector{T}, x::AbstractVector{T}, q95::T) where {T<:Real}
  rollmax5!(x5, x)
  rx5day = maximum(x5)

  r20mm = sum(x .>= 20.0) # days
  r95ptot = sum(x[x.>q95])
  prcptot = sum(x[x.>=1.0])
  cdd = CDD(x)
  [cdd, r20mm, rx5day, r95ptot, prcptot]
end


function index_P(x::AbstractVector{T}, q95::T) where {T<:Real}
  x5 = rollmax(x, 5)
  rx5day = maximum(x5)

  r20mm = sum(x .>= 20.0) # days
  r95ptot = sum(x[x.>q95])
  prcptot = sum(x[x.>=1.0])
  cdd = CDD(x)
  [cdd, r20mm, rx5day, r95ptot, prcptot]
end
