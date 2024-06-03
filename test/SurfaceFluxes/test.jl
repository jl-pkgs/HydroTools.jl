@time using Roots
f(x) = exp(x) - x^4

## bracketing
fzero(f, 8, 9)          # 8.613169456441398
fzero(f, -10, 0)      # -0.8155534188089606
fzeros(f, -10, 10)

# Brent只能找到一个
fzero(f, -10, 10, Roots.Brent())          # 8.6131694564414

@run root_brent(f, lb=8, ub=9, tol=0.001)
