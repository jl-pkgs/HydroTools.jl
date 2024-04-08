export selectEventMaxDuration, extract_eventInfo,
  detect_events


# import StatsBase: rle, inverse_rle
which_isna(x) = findall(x .== nothing)
which_notna(x) = findall(x .!= nothing)

function rm_empty(x)
  inds = which_notna(x)
  x[inds]
end

function selectEventMaxDuration(events::AbstractVector)
  # @show events
  if length(events) > 0
    _, pos = findmax(map(x -> x.len, events))
    events[pos]
  end
end


function detect_events(lgl::AbstractVector{Bool};
  len_min=1,
  index=nothing, only_max=false, ignored...)

  I = findall(lgl)
  flag = cumsum([1; diff(I) .!= 1])
  n = flag[end]

  events = map(k -> begin
      _ind = I[flag.==k]
      i_beg = _ind[1]
      i_end = _ind[end]

      len = i_end - i_beg + 1

      if index !== nothing
        i_beg = index[i_beg]
        i_end = index[i_end]
      end

      if len >= len_min
        (; i_beg, i_end, len)
      end # else nothing
    end, 1:n) |> rm_empty
  # events = events[which_notna(events)]
  only_max && (events = selectEventMaxDuration(events))
  events
  # use_df ? DataFrame(events) : events
end


"""
    detect_events(lgl::AbstractVector{Bool}; len_min=1, index=nothing, only_max=false, ignored...)
    detect_events(y::AbstractVector{<:Real}, lgl::BitVector;
        len_min=1,
        len2peak=0,
        goal=1,
        index=nothing,
        only_max=false,
        ignored...
    )

Detects events in a signal based on a logical vector.

# Arguments
- `y::AbstractVector{<:Real}`: The signal to detect events in.
- `lgl::BitVector`: A logical vector indicating where events occur in the signal.
- `len_min=1`: (optional) The minimum length of an event.
- `len2peak=0`: (optional) The minimum length from the peak to the beginning or end of an event.
- `goal=1`: (optional), `-1` or `1`. If `1`, find the maximum value in the event; If `-1`, minimum value used.
- `index=nothing`: (optional) The index of the signal. If provided, the returned indices will be in the index space.
- `only_max=false`: (optional) If `true`, only events with the maximum duration will be returned.
- `ignored...`: (optional) Ignored arguments.

# Returns
An array of named tuples, where each tuple represents an event and has the following fields:
- `i_beg::Integer`: The index of the beginning of the event.
- `i_peak::Integer`: The index of the peak of the event.
- `i_end::Integer`: The index of the end of the event.
- `len::Integer`: The length of the event.
- `len_left::Integer`: The length from the beginning of the event to the peak.
- `len_right::Integer`: The length from the peak to the end of the event.
- `peak::Real`: The value of the peak of the event.

If `only_max` is `true`, only the event with the maximum duration will be returned.

# Examples
```julia
julia> y = [0, 1, 2, 3, 0, 2, 1, 0];
julia> lgl = y .> 0;
julia> detect_events(lgl)
2-element Vector{NamedTuple{(:i_beg, :i_end, :len), Tuple{Int64, Int64, Int64}}}:
 (i_beg = 2, i_end = 4, len = 3)
 (i_beg = 6, i_end = 7, len = 2)
julia> detect_events(y, lgl)
2-element Vector{NamedTuple{(:i_beg, :i_peak, :i_end, :len, :len_left, :len_right, :peak), NTuple{7, Int64}}}:
 (i_beg = 2, i_peak = 4, i_end = 4, len = 3, len_left = 2, len_right = 0, peak = 3)
 (i_beg = 6, i_peak = 6, i_end = 7, len = 2, len_left = 0, len_right = 1, peak = 2)
```
"""
function detect_events(y::AbstractVector{<:Real}, lgl::AbstractVector{Bool};
  len_min=1, len2peak=0, goal=1,
  index=nothing, only_max=false, ignored...)

  I = findall(lgl)
  flag = cumsum([1; diff(I) .!= 1])
  n = flag[end]

  events = map(k -> begin
      _ind = I[flag.==k]
      i_beg = _ind[1]
      i_end = _ind[end]

      event = extract_eventInfo(y, i_beg, i_end; index, len_min, len2peak, goal)
    end, 1:n) |> rm_empty

  only_max && (events = selectEventMaxDuration(events))
  events
  # use_df ? DataFrame(events) : events
end



"""
    extract_eventInfo(y::AbstractVector{<:Real}, i_beg::Integer, i_end::Integer;
        index=nothing, goal=1, len_min=1, len2peak=0, ignored...)

Extracts information about an event in a signal.

# Examples
```julia
julia> y = [0, 1, 2, 3, 2, 1, 0];
julia> extract_eventInfo(y, 2, 6)
(i_beg = 2, i_peak = 4, i_end = 5, len = 4, len_left = 2, len_right = 1, peak = 3)
```
"""
function extract_eventInfo(y::AbstractVector{<:Real}, i_beg::Integer, i_end::Integer;
  index=nothing,
  goal=1,
  len_min=1, len2peak=0, ignored...)

  inds = i_beg:i_end
  yi = @view y[inds]

  if goal == 1
    peak, _i_peak = findmax(yi)
  elseif goal == -1
    peak, _i_peak = findmin(yi)
  end

  i_peak = inds[_i_peak]

  len = i_end - i_beg + 1
  len_right = i_end - i_peak
  len_left = i_peak - i_beg

  if index !== nothing
    i_beg = index[i_beg]
    i_end = index[i_end]
    i_peak = index[i_peak]
  end

  if len >= len_min && min(len_right, len_left) >= len2peak
    (; i_beg, i_peak, i_end, len, len_left, len_right, peak)
  end # else nothing
end

# function zip_continue(x::AbstractVector{Bool})
#   flag = cumsum([1; diff(x) .!= 0])
#   grps = unique(flag)
#   n = grps[end]

#   res = []
#   for i = 1:n
#     inds = findall(flag .== grps[i])
#     index = [inds[1], inds[end]]
#     push!(res, (; grp=i, index, value=x[index[1]]))
#   end
#   res
# end
# export zip_continue;
