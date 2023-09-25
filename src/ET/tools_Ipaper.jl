using Statistics


# main script of moving average
function movmean2(y::AbstractVector{T}, win_left::Integer, win_right::Integer=win_right) where {T<:Real}
  n = length(y)
  z = zeros(Float64, n)

  @inbounds for i in 1:n
    i_beg = max(i - win_left, 1)
    i_end = min(i + win_right, n)

    n_i = 0        # number
    sum = T(0.0)   # sum of values in window

    for j in i_beg:i_end
      n_i += 1
      sum += y[j]
    end
    z[i] = sum / n_i
  end
  z
end


function nanmean2(x, y)
  if isnan(x)
    y
  elseif isnan(y)
    x
  else
    (x + y) / 2.0
  end
end


function replace_miss(x::AbstractVector, miss=NaN)
  x[ismissing.(x)] .= miss
  x
end

function getDataType(x)
  type = eltype(x)
  typeof(type) == Union ? type.b : type
end


function struct2vec(x)
  keys = fieldnames(typeof(x))
  vals = [getfield(x, key) for key in keys]
  vals
end

function struct2tuple(x)
  keys = fieldnames(typeof(x))
  vals = [getfield(x, key) for key in keys]
  (; zip(keys, vals)...)
end


export struct2vec, struct2tuple
export movmean2, nanmean2, getDataType, replace_miss
