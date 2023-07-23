"""
    cceua(fn, s, sf, bl, bu, icall)

Generate a new point in a simplex

# Example
```julia
snew, fnew, icall = cceua(fn, s, sf, bl, bu, icall)
```
"""
function cceua(fn, s, sf, bl, bu, icall)
    #  This is the subroutine for generating a new point in a simplex
    #
    #   s[.,.] = the sorted simplex in order of increasing function values
    #   s[.] = function values in increasing order
    #
    # LIST OF LOCAL VARIABLES
    #   sb[.] = the best point of the simplex
    #   sw[.] = the worst point of the simplex
    #   w2[.] = the second worst point of the simplex
    #   fw = function value of the worst point
    #   ce[.] = the centroid of the simplex excluding wo
    #   snew[.] = new point generated from the simplex
    #   iviol = flag indicating if constraints are violated
    #         = 1 ; yes
    #         = 0 ; no

    nps, nopt = size(s)
    n = nps
    # m = nopt
    alpha = 1.0
    beta = 0.5

    # Assign the best & worst points:
    # sb = s[1, :]
    # fb = sf[1]
    sw = s[n, :]
    fw = sf[n]

    # Compute the centroid of the simplex excluding the worst point:
    ce = mean(s[1:n-1, :], dims = 1)[:]

    # Attempt a reflection point
    snew = ce .+ alpha * (ce .- sw)

    # Check if is outside the bounds:
    ibound = 0
    
    s1 = snew - bl
    idx = findall(s1 .< 0)
    !isempty(idx) && (ibound = 1)
    
    s1 = bu - snew
    idx = findall(s1 .< 0)
    !isempty(idx) && (ibound = 2)
    
    if ibound >= 1
        # snew = bl + rand(nopt) .* (bu - bl)
        snew = bl + mrand(nopt) .* (bu - bl)
    end

    fnew = fn(snew)
    icall = icall + 1

    # Reflection failed; now attempt a contraction point:
    if fnew .> fw
        snew = sw .+ beta * (ce .- sw)
        fnew = fn(snew)
        icall = icall + 1

        # Both reflection & contraction have failed; attempt a random point
        if fnew .> fw
            # snew = bl + rand(nopt) .* (bu - bl)
            snew = bl + mrand(nopt) .* (bu - bl)
            fnew = fn(snew)
            icall = icall + 1
        end
    end
    snew, fnew, icall
end
