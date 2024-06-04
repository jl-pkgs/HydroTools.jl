"""
Solve for the root of a function using the secant and Brent's methods given
initial estimates xa and xb. The root is updated until its accuracy is tol. func
is the name of the function to solve. The variable root is returned as the root
of the function. The function being evaluated has the definition statement:

function [fx] = func (physcon, forcvar, surfvar, x)

The function func is evaluated at x and the returned value is fx. It uses
variables in the physcon, forcvar, surfvar, and fluxvar structures. These are
passed in as input arguments. It also calculates values for variables in the
fluxvar structure so this must be returned in the function call as an output
argument. The Julia function `func` evaluates func.

# INPUTS
- `args...`: other arguments to be passed to `func`
- `tol`: accuracy tolerance for x
"""
function root_hybrid(func, args...; lb::Real, ub::Real, tol=0.01, itmax = 40, kw...)
  # --- Evaluate func at `xa` and `xb` and see if this is the root
  x0 = lb
  f0 = func(x0, args...; kw...)
  (f0 == 0) && return x0
  
  x1 = ub
  f1 = func(x1, args...; kw...)
  (f1 == 0) && return x1
  
  # --- Order initial root estimates correctly
  if (f1 < f0)
    minx, minf = x1, f1
  else
    minx, minf = x0, f0
  end

  # --- Iterative root calculation. Use the secant method, with Brent's method as a backup
  for iter = 1:itmax
    dx = -f1 * (x1 - x0) / (f1 - f0)
    x = x1 + dx

    # Check if x is the root. If so, exit the iteration
    if (abs(dx) < tol)
      x0 = x
      break
    end

    # Evaluate the function at x
    x0, f0 = x1, f1
    x1 = x
    f1 = func(x1, args...; kw...)
    if (f1 < minf)
      minx, minf = x1, f1
    end

    # If a root zone is found, use Brent's method for a robust backup strategy
    # and exit the iteration
    if (f1 * f0 < 0)
      x = root_brent(func, args...; lb=x0, ub=x1, tol, kw...)
      x0 = x
      break
    end

    # In case of failing to converge within itmax iterations stop at the minimum function
    if (iter == itmax)
      f1 = func(minx, args...; kw...)
      x0 = minx
    end
  end
  return x0
end
