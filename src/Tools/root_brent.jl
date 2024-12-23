"""
Use Brent's method to find the root of a function, which is known to exist
between xa and xb. The root is updated until its accuracy is tol. func is the
name of the function to solve. The variable root is returned as the root of
the function. The function being evaluated has the definition statement: 

`function [fx] = func (physcon, forcvar, surfvar, x)`

The function func is evaluated at x and the returned value is fx. It uses
variables in the physcon, forcvar, surfvar, and fluxvar structures. These are
passed in as input arguments. It also calculates values for variables in the
fluxvar structure so this must be returned in the function call as an output
argument. The Julia function `func` evaluates func.

# INPUTS
- `args...`: other arguments to be passed to `func`
"""
function root_brent(func, args...; lb, ub, tol=0.01, eps=1e-8, itmax = 50, kw...)
  # --- Evaluate func at xa and xb and make sure the root is bracketed
  
  a, b = lb, ub
  fa = func(a, args...; kw...)
  fb = func(b, args...; kw...)

  if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0))
    error("root_brent error: root must be bracketed")
  end
  
  # Relative error tolerance
  c, fc = b, fb
  d = 0.0
  e = 0.0
  # --- Iterative root calculation
  for iter = 1:itmax
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
      c, fc = a, fa
      d = b - a
      e = d
    end
    if (abs(fc) < abs(fb))
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
    end
    tol1 = 2 * eps * abs(b) + 0.5 * tol
    xm = 0.5 * (c - b)

    # Check to end iteration
    (abs(xm) <= tol1 || fb == 0) && break

    if (abs(e) >= tol1 && abs(fa) > abs(fb))
      s = fb / fa
      if (a == c)
        p = 2 * xm * s
        q = 1 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1))
        q = (q - 1) * (r - 1) * (s - 1)
      end
      if (p > 0)
        q = -q
      end
      p = abs(p)
      min_val = min(3 * xm * q - abs(tol1 * q), abs(e * q))
      if (2 * p < min_val)
        e = d
        d = p / q
      else
        d = xm
        e = d
      end
    else
      d = xm
      e = d
    end
    a, fa = b, fb
    if (abs(d) > tol1)
      b = b + d
    else
      b = xm >= 0 ? b + abs(tol1) : b - abs(tol1)
    end

    fb = func(b, args...; kw...)
    fb == 0 && break # Check to end iteration

    # Check to see if failed to converge
    (iter == itmax) && error("root_brent error: Maximum number of iterations exceeded")
  end
  return b
end
