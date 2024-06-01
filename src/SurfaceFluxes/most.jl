"""
Use Monin-Obukhov similarity theory to obtain the Obukhov length (obu).

This is the function to solve for the Obukhov length. For the current
estimate of the Obukhov length (x), calculate u*, T*, and q* and then
the new length (obu). The function value is the change in Obukhov length:
fx = x - obu.
"""
function most(x, physcon, forcvar, fluxvar)
  # physcon, forcvar, _, fluxvar = flatten(varargin)
  # Prevent near-zero values of the Obukhov length
  x = abs(x) <= 0.1 ? 0.1 : x

  # Calculate z-d at the reference height, because this is used many times
  z_minus_d = forcvar.zref - fluxvar.disp # z - d

  # Evaluate psi for momentum at the reference height (zref-disp) and surface (z0m)
  psi_m_zref = psi_m_monin_obukhov(z_minus_d / x) # x = L_{MO}
  psi_m_z0m = psi_m_monin_obukhov(fluxvar.z0m / x)
  psim = -psi_m_zref + psi_m_z0m

  # Evaluate psi for scalars at the reference height (zref-disp) and surface (z0c)
  psi_c_zref = psi_c_monin_obukhov(z_minus_d / x)
  psi_c_z0c = psi_c_monin_obukhov(fluxvar.z0c / x)
  psic = -psi_c_zref + psi_c_z0c

  # Calculate u* (m/s), T* (K), q* (mol/mol), and Tv* (K)
  zlog_m = log(z_minus_d / fluxvar.z0m)
  zlog_c = log(z_minus_d / fluxvar.z0c)

  fluxvar.ustar = forcvar.uref * physcon.vkc / (zlog_m + psim)
  fluxvar.tstar = (forcvar.thref - fluxvar.tsrf) * physcon.vkc / (zlog_c + psic)
  fluxvar.qstar = (forcvar.eref - fluxvar.esrf) / forcvar.pref * physcon.vkc / (zlog_c + psic)
  tvstar = fluxvar.tstar + 0.61 * forcvar.thref * fluxvar.qstar * (physcon.mmh2o / forcvar.mmair)

  # Calculate Obukhov length (m)
  fluxvar.obu = fluxvar.ustar^2 * forcvar.thvref / (physcon.vkc * physcon.grav * tvstar)

  # Calculate change in obu
  fx = x - fluxvar.obu
  return fx
end
