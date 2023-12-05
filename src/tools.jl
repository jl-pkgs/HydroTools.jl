"""
    f_tridiagonal_solver(a, b, c, d, n::Integer)

Tridiagonal solver

# Description

Converted into a R code from the original code of Gordon Bonan: Bonan, G.
(2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge:
Cambridge University Press. doi:10.1017/9781107339217

Solve for U given the set of equations R * U = D, where U is a vector
of length N, D is a vector of length N, and R is an N x N tridiagonal
matrix defined by the vectors A, B, C each of length N. A(1) and
C(N) are undefined and are not referenced.

    |B(1) C(1) ...  ...  ...                     |
    |A(2) B(2) C(2) ...  ...                     |
R = |     A(3) B(3) C(3) ...                     |
    |                    ... A(N-1) B(N-1) C(N-1)|
    |                    ... ...    A(N)   B(N)  |

The system of equations is written as:

  A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i

for i = 1 to N. The solution is found by rewriting the
equations so that:

  U_i = F_i - E_i * U_i+1

# Return
- `Solution: U`
"""
function f_tridiagonal_solver(a, b, c, d, n::Integer)
  e = fill(NaN, n - 1)
  e[1] = c[1] / b[1]

  for i in 2:(n-1)
    e[i] = c[i] / (b[i] - a[i] * e[i-1])
  end

  f = fill(NaN, n)
  f[1] = d[1] / b[1]

  for i in 2:(n)
    f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
  end

  u = fill(NaN, n)
  u[n] = f[n]

  for i in n-1:-1:1
    u[i] = f[i] - e[i] * u[i+1]
  end
  return (u)
end

export f_tridiagonal_solver
