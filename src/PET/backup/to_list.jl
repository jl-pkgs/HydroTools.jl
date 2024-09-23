import LabelledArrays
import LabelledArrays: LVector, LArray, symnames

const LA = LabelledArrays

LA.LVector(keys::Vector{Symbol}, values) = LVector(; zip(keys, values)...)
LA.LVector(keys::Vector{<:AbstractString}, values) = LVector(; zip(Symbol.(keys), values)...)
LA.LVector(keys::Tuple, values) = LVector(; zip(keys, values)...)

LA.LVector(keys::Vector{Symbol}) = LVector(keys, zeros(length(keys)))
LA.LVector(keys::Tuple) = LVector(keys, zeros(length(keys)))


list = LVector;
Base.names(x::LArray) = symnames(typeof(x))


function add(x::LArray, y::LArray)
  list([keys(x)..., keys(y)...],
    [values(x)..., values(y)...],)
end

export list, add
