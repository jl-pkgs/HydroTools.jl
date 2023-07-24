export GroupFlag, CDD, CWD, index_P!, index_P


function GroupFlag(inds::Vector{Int64})
  cumsum([1; diff(inds) .!= 1])
end

# using RollingFunctions
# using StatsBase

include("index_P.jl")
# R20, Rx5day, R95pTOT, CDD, prcptot

