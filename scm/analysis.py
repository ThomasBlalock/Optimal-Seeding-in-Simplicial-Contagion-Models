import numpy as np


def mf_rho(lam, lam_d):
    """Mean-field analytical steady state of the SCM.

    Returns the real positive root of:
        rho = (lam_d - lam +/- sqrt((lam - lam_d)^2 - 4*lam_d*(1 - lam))) / (2*lam_d)
    or 0 when the discriminant is negative (below threshold).
    """
    if lam_d == 0:
        return 1 - 1 / lam

    discriminant = (lam - lam_d) ** 2 - 4 * lam_d * (1 - lam)
    sqrt_term = np.sqrt(discriminant.astype(complex))
    rho = (lam_d - lam + sqrt_term) / (2 * lam_d)
    rho = np.where(np.imag(rho) != 0, 0, rho)

    return rho.real


def steady_state_rho(rho_history, t_avg=100):
    """Compute steady-state infection density from a simulation history."""
    if rho_history[-1] == 0.0:
        return 0.0
    return np.mean(rho_history[-t_avg:])
