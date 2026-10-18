"""Finite-volume solver for steady fully developed flow of a generalised Newtonian fluid.

Governing equation (half channel 0 <= y <= h, by symmetry):

    d/dy [ mu(|du/dy|) du/dy ] = -G ,     du/dy(0) = 0 ,   u(h) = 0

Discretisation: uniform grid, conservative second-order scheme. Face viscosity is
evaluated from the face shear rate (u_{i+1}-u_i)/dy. The non-linearity is resolved by
Picard (successive-substitution) iteration with under-relaxation; each linear
tridiagonal system is solved directly (banded LU).
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.linalg import solve_banded

_trapz = getattr(np, "trapezoid", None) or np.trapz


@dataclass
class ChannelResult:
    y: np.ndarray
    u: np.ndarray
    G: float
    h: float
    iterations: int
    converged: bool
    residual_history: list = field(default_factory=list)

    @property
    def bulk_velocity(self):
        return _trapz(self.u, self.y) / self.h

    def shear_rate(self):
        """|du/dy| at nodes: 2nd-order central inside, 2nd-order one-sided at the wall."""
        u = self.u
        dy = self.y[1] - self.y[0]
        g = np.empty_like(u)
        g[1:-1] = (u[2:] - u[:-2]) / (2 * dy)
        g[0] = 0.0  # symmetry
        g[-1] = (3 * u[-1] - 4 * u[-2] + u[-3]) / (2 * dy)
        return np.abs(g)

    def wall_shear_rate(self):
        return float(self.shear_rate()[-1])

    def wall_shear_stress(self, viscosity):
        gw = self.wall_shear_rate()
        return float(viscosity(np.array([gw]))[0] * gw)


def solve_channel(N, G, viscosity, h=1.0, tol=1e-10, max_iter=20000, relax=0.5):
    """Solve for u(y) on N intervals of the half channel [0, h].

    viscosity : callable mu(gamma_dot) -> apparent viscosity (vectorised).
    """
    y = np.linspace(0.0, h, N + 1)
    dy = y[1] - y[0]
    M = N  # unknowns u_0..u_{N-1}; u_N = 0 (no slip)

    mu_ref = float(viscosity(np.array([1.0]))[0])
    u = G / (2 * mu_ref) * (h**2 - y**2)  # Newtonian-like initial guess
    hist, converged = [], False
    it = 0
    for it in range(1, max_iter + 1):
        mf = viscosity(np.abs(np.diff(u)) / dy)       # viscosity at the N faces
        main = np.empty(M)
        main[0] = -2.0 * mf[0]                        # symmetry node
        main[1:] = -(mf[1:M] + mf[0:M - 1])
        ab = np.zeros((3, M))
        ab[1, :] = main
        ab[0, 1] = 2.0 * mf[0]                        # row 0 upper coefficient
        ab[0, 2:] = mf[1:M - 1]                       # rows i>=1 upper coefficients
        ab[2, :-1] = mf[0:M - 1]                      # rows i>=1 lower coefficients
        rhs = np.full(M, -G * dy**2)
        u_new = np.append(solve_banded((1, 1), ab, rhs), 0.0)
        u_next = (1 - relax) * u + relax * u_new
        res = float(np.max(np.abs(u_next - u)) / np.max(np.abs(u_next)))
        hist.append(res)
        u = u_next
        if res < tol:
            converged = True
            break
    return ChannelResult(y=y, u=u, G=G, h=h, iterations=it, converged=converged,
                         residual_history=hist)
