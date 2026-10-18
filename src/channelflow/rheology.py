"""Constitutive models: apparent viscosity mu(gamma_dot) for generalised Newtonian fluids.

All functions take the (non-negative) shear-rate magnitude ``gd`` and return the
apparent viscosity, so that shear stress tau = mu(gd) * gd.
"""
import numpy as np


def newtonian_viscosity(gd, mu=1.0):
    return mu * np.ones_like(np.asarray(gd, dtype=float))


def power_law_viscosity(gd, K=1.0, n=1.0, eps=1e-8):
    """Ostwald-de Waele: mu = K * gd**(n-1).

    For n < 1 the viscosity diverges as gd -> 0 (on the channel centreline).
    A tiny regularisation ``eps`` is added: mu = K * (gd^2 + eps^2)^((n-1)/2),
    which leaves the solution unchanged to well below discretisation error
    while keeping the discrete system finite.
    """
    gd = np.asarray(gd, dtype=float)
    return K * (gd**2 + eps**2) ** ((n - 1.0) / 2.0)


def carreau_viscosity(gd, mu0=1.0, mu_inf=0.01, lam=1.0, n=0.5):
    """Carreau model: mu = mu_inf + (mu0 - mu_inf) * [1 + (lam*gd)^2]^((n-1)/2).

    Bounded at both ends (Newtonian plateaus), so it removes the unphysical
    infinite zero-shear viscosity of the power law.
    """
    gd = np.asarray(gd, dtype=float)
    return mu_inf + (mu0 - mu_inf) * (1.0 + (lam * gd) ** 2) ** ((n - 1.0) / 2.0)
