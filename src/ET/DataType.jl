using Parameters


@with_kw mutable struct interm_PML{T}
  ET::T = T(0)
  GPP::T = T(0)
  Ec::T = T(0)
  Ecr::T = T(0)
  Eca::T = T(0)

  Ei::T = T(0)
  Pi::T = T(0)
  Es_eq::T = T(0)
  Eeq::T = T(0)
  ET_water::T = T(0)

  Ga::T = T(0)
  Gc_w::T = T(0)

  fval_soil::T = T(0)
  Es::T = T(0)
end

@with_kw mutable struct output_PML{T}
  n::Integer
  ET::Vector{T} = zeros(T, n)
  GPP::Vector{T} = zeros(T, n)
  Ec::Vector{T} = zeros(T, n)
  Ecr::Vector{T} = zeros(T, n)
  Eca::Vector{T} = zeros(T, n)

  Ei::Vector{T} = zeros(T, n)
  Pi::Vector{T} = zeros(T, n)
  Es_eq::Vector{T} = zeros(T, n)
  Eeq::Vector{T} = zeros(T, n)
  ET_water::Vector{T} = zeros(T, n)

  Ga::Vector{T} = zeros(T, n)
  Gc_w::Vector{T} = zeros(T, n)

  fval_soil::Vector{T} = zeros(T, n)
  Es::Vector{T} = zeros(T, n)
end
# output_PML(;n::Integer) = output_PML{Float64}(;n=n)


Base.getindex(x::interm_PML, key::Union{Symbol,String}) = getfield(x, Symbol(key))

function Base.setindex!(res::output_PML, r::Union{NamedTuple,interm_PML}, t, fields)
  # fields = fieldnames(output_PML)[2:end]
  for i = eachindex(fields)
    field = fields[i]
    x = getfield(res, field)
    x[t] = r[field]
  end
end

# function setindex!(res::output_PML, r::Union{NTuple,Vector}, t, fields)
function Base.setindex!(res::output_PML, r::NTuple, t, fields)
  for i = eachindex(fields)
    field = fields[i]
    x = getfield(res, field)
    x[t] = r[i]
  end
end


## DATATYPE CONVERSION ---------------------------------------------------------

function to_mat(res::output_PML{T}) where {T<:Real}
  names = fieldnames(output_PML)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data, names
end


export interm_PML, output_PML
export to_mat;

