"""Closed-form solutions for steady, fully developed, pressure-driven flow between
parallel plates at y = +/- h with constant pressure gradient G = -dp/dx > 0.

Momentum balance (fully developed):  d(tau)/dy = -G  ->  tau(y) = G*y  (symmetry).
"""
import numpy as np


def poiseuille_velocity(y, G=1.0, mu=1.0, h=1.0):
    """Newtonian plane Poiseuille flow: u = G/(2 mu) (h^2 - y^2)."""
    y = np.asarray(y, dtype=float)
    return G / (2.0 * mu) * (h**2 - y**2)


def powerlaw_velocity(y, G=1.0, K=1.0, n=1.0, h=1.0):
    """Power-law fluid, tau = K |du/dy|^(n-1) du/dy:

        u(y) = n/(n+1) (G/K)^(1/n) [ h^((n+1)/n) - |y|^((n+1)/n) ]
    """
    y = np.abs(np.asarray(y, dtype=float))
    p = (n + 1.0) / n
    return n / (n + 1.0) * (G / K) ** (1.0 / n) * (h**p - y**p)


def bulk_velocity_powerlaw(G=1.0, K=1.0, n=1.0, h=1.0):
    """Mean velocity Ub = (1/h) int_0^h u dy = n/(2n+1) (G/K)^(1/n) h^((n+1)/n)."""
    return n / (2.0 * n + 1.0) * (G / K) ** (1.0 / n) * h ** ((n + 1.0) / n)


def pressure_gradient_for_bulk(Ub, K=1.0, n=1.0, h=1.0):
    """Invert Ub(G): the pressure gradient needed to drive a given bulk velocity."""
    return K * (Ub * (2.0 * n + 1.0) / (n * h ** ((n + 1.0) / n))) ** n


def wall_shear_rate_powerlaw(Ub, n=1.0, h=1.0):
    """|du/dy| at the wall = (2n+1)/n * Ub/h   (n=1 gives the Newtonian 3 Ub/h)."""
    return (2.0 * n + 1.0) / n * Ub / h
