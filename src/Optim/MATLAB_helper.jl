disp = println;
num2str = string;

function SORT(x::Vector{<:Real})
  idx = sortperm(x)
  x[idx], idx
end

function MAX(x, dims=1)
  maximum(x, dims=dims)[:]
end

function MIN(x, dims=1)
  minimum(x, dims=dims)[:]
end

function MEAN(x, dims)
  mean(x, dims=dims)[:]
end
MEAN(x) = MEAN(x, 1)


## RANDOM ----------------------------------------------------------------------
set_seed(seed=1) = Random.seed!(seed)
mrand() = rand()
mrand(n) = rand(n)

# using MATLAB
# set_seed(seed=1) = mat"rand('seed', $seed);"
# mrand() = mat"rand();"
# mrand(n) = mat"rand($n, 1);"

export SORT, MEAN, MAX, MIN
