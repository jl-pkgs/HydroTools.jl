T0 = 273.15

K0 = 273.15

L = 2.5e6

Es_T0 = 6.11

atm = 101.325

#' Molecular weight
#' 
#' @description
#' `Mw`: Molecular weight of dry air \eqn{M_d = 28.9634g/mol}
#' `Md`: [Molecular weight](https://en.wikipedia.org/wiki/Molar_mass) of water vapor \eqn{M_w = 18.01528g/mol}
#' `epsilon`: the ratio of `Mw` to `Md`
Mw = 18.01528

#' @rdname Mw
Md = 28.9634

#' @rdname Mw
epsilon = Mw / Md

#' Specific gas constants
#' 
#' @description
#' - `R`: [gas constant](https://en.wikipedia.org/wiki/Gas_constant), J/(mol K)
#' - `Rw`: [Specific gas constant](https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant)
#' of water vapor \eqn{R_w = \frac{1000R}{M_w} = 461.52J/(kg K)}.
#' - `Rd`: Specific gas constant of dry air (J/(kg K)).
#' @export
R = 8.3144621 # J/(mol K)

#' @rdname R
#' @export
Rw = R / Mw * 1000. # J/(kg K)

#' @rdname R
#' @export
Rd = R / Md * 1000.

Cp = 1.013 * 1e-3 # MJ kg-1 degC-1

