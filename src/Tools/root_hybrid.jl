"""
Solve for the root of a function using the secant and Brent's methods given
initial estimates xa and xb. The root is updated until its accuracy is tol. func
is the name of the function to solve. The variable root is returned as the root
of the function. The function being evaluated has the definition statement:

function [fluxvar, fx] = func (physcon, forcvar, surfvar, fluxvar, x)

The function func is evaluated at x and the returned value is fx. It uses
variables in the physcon, forcvar, surfvar, and fluxvar structures. These are
passed in as input arguments. It also calculates values for variables in the
fluxvar structure so this must be returned in the function call as an output
argument. The Julia function `func` evaluates func.

# INPUTS
- `args...`: other arguments to be passed to `func`
"""
function root_hybrid(func, xa, xb, tol, args...; kw...)
  # --- Evaluate func at xa and see if this is the root
  x0 = xa
  fluxvar, f0 = func(x0, args...; kw...)
  if (f0 == 0)
    return fluxvar, x0
  end

  # --- Evaluate func at xb and see if this is the root
  x1 = xb
  fluxvar, f1 = func(x1, args...; kw...)
  if (f1 == 0)
    return fluxvar, x1
  end

  # --- Order initial root estimates correctly
  if (f1 < f0)
    minx, minf = x1, f1
  else
    minx, minf = x0, f0
  end

  # --- Iterative root calculation. Use the secant method, with Brent's method as a backup
  itmax = 40
  for iter = 1:itmax
    dx = -f1 * (x1 - x0) / (f1 - f0)
    x = x1 + dx

    # Check if x is the root. If so, exit the iteration
    if (abs(dx) < tol)
      x0 = x
      break
    end

    # Evaluate the function at x
    x0 = x1
    f0 = f1
    x1 = x
    fluxvar, f1 = func(x1, args...; kw...)
    if (f1 < minf)
      minx, minf = x1, f1
    end

    # If a root zone is found, use Brent's method for a robust backup strategy
    # and exit the iteration
    if (f1 * f0 < 0)
      fluxvar, x = root_brent(func, x0, x1, tol, args...; kw...)
      x0 = x
      break
    end

    # In case of failing to converge within itmax iterations stop at the minimum function
    if (iter == itmax)
      fluxvar, f1 = func(minx, args...; kw...)
      x0 = minx
    end
  end

  return fluxvar, x0
end
