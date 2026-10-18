"""Error norms and grid-convergence tools."""
import numpy as np


def l2_error(u_num, u_exact):
    """Relative discrete L2 error."""
    return float(np.linalg.norm(u_num - u_exact) / np.linalg.norm(u_exact))


def linf_error(u_num, u_exact):
    """Relative maximum-norm error."""
    return float(np.max(np.abs(u_num - u_exact)) / np.max(np.abs(u_exact)))


def observed_order(e_coarse, e_fine, r=2.0):
    """Observed order of accuracy p = ln(e_c/e_f)/ln(r) for refinement ratio r."""
    if e_fine <= 0 or e_coarse <= 0:
        return float("nan")
    return float(np.log(e_coarse / e_fine) / np.log(r))


def gci(f_fine, f_medium, f_coarse, r=2.0, Fs=1.25):
    """Grid Convergence Index (Roache) from three systematically refined grids.

    Returns (p, f_extrapolated, GCI_fine) with GCI as a fraction of f_fine.
    """
    e32 = f_coarse - f_medium
    e21 = f_medium - f_fine
    if e21 == 0 or e32 == 0:
        return float("nan"), f_fine, 0.0
    p = np.log(abs(e32 / e21)) / np.log(r)
    f_ext = f_fine + (f_fine - f_medium) / (r**p - 1.0)
    ea = abs((f_fine - f_medium) / f_fine)
    return float(p), float(f_ext), float(Fs * ea / (r**p - 1.0))
