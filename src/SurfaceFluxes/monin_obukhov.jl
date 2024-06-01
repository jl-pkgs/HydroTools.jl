# ϕ: phi
# ψ: psi
function phi_m_monin_obukhov(x)
  # Evaluate the Monin-Obukhov phi function for momentum at x
  # Eq. 6.37
  if x < 0
    phi_m = (1 - 16 * x)^(-0.25)
  else
    phi_m = 1 + 5 * x
  end
  return phi_m
end

function phi_c_monin_obukhov(x)
  # Evaluate the Monin-Obukhov phi function for scalars at x
  # Eq. 6.38
  if x < 0
    phi_c = (1 - 16 * x)^(-0.5)
  else
    phi_c = 1 + 5 * x
  end
  return phi_c
end

function psi_c_monin_obukhov(x)
  # Evaluate the Monin-Obukhov ψ function for scalars at x
  # Eq. 6.47
  if x < 0
    y = (1 - 16 * x)^0.25
    psi_c = 2 * log((1 + y^2) / 2)
  else
    psi_c = -5 * x
  end
  return psi_c
end

function psi_m_monin_obukhov(x)
  # Evaluate the Monin-Obukhov ψ function for momentum at x
  # Eq. 6.46
  if x < 0
    y = (1 - 16 * x)^0.25
    psi_m = 2 * log((1 + y) / 2) + log((1 + y^2) / 2) - 2 * atan(y) + pi / 2
  else
    psi_m = -5 * x
  end
  return psi_m
end
