# Canopy conductance (mol/m2/s) - use a weighted average of sunlit and shaded leaves
gc_min = 0.05                           # Minimum conductance (mol/m2/s)
gc_max = 0.2                            # Maximum conductance (mol/m2/s)

"""
# Fields
$(TYPEDFIELDS)

# Example
```julia
coszen = cal_coszen(doy, hour, lat)
```
"""
@with_kw mutable struct Canopy{T<:AbstractFloat}
  "sin of solar zenith angle"
  coszen::T = T(NaN)
  "Leaf area index [m2 m-2]"
  LAI::T = 2.0
  "Light extinction coefficient"
  kQ::T = 0.5 / max(coszen, 0.0001)
  "Sunlit fraction of canopy"
  fsun::T = (1 - exp(-kQ * LAI)) / (kQ * LAI) # 大叶模型
  "canopy conductance (mol/m2/s)"
  gc::T = coszen > 0.0 ? (fsun * gc_max + (1 - fsun) * gc_min) * LAI : gc_min * LAI # 冠层整体的导度
end

function update_canopy!(can::Canopy{T}, LAI::T, coszen::T) where {T<:Real}
  (; kQ) = can
  fsun = (1 - exp(-kQ * LAI)) / (kQ * LAI)
  gc = coszen > 0.0 ? (fsun * gc_max + (1 - fsun) * gc_min) * LAI : gc_min * LAI
  @pack! can = LAI, coszen, fsun, gc
  can
end

function update_canopy!(can::Canopy{T}, coszen::T) where {T<:Real}
  (; fsun, LAI) = can
  gc = coszen > 0.0 ? (fsun * gc_max + (1 - fsun) * gc_min) * LAI : gc_min * LAI
  @pack! can = coszen, gc
  can
end


export Canopy, update_canopy!
