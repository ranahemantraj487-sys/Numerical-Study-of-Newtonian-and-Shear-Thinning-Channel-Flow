import numpy as np
import pytest
from channelflow import *


def test_newtonian_matches_poiseuille():
    r = solve_channel(40, 2.0, lambda g: newtonian_viscosity(g, mu=0.5))
    ex = poiseuille_velocity(r.y, G=2.0, mu=0.5)
    assert r.converged
    assert l2_error(r.u, ex) < 1e-10          # parabola is exact for 2nd-order FV


@pytest.mark.parametrize("n", [0.9, 0.7, 0.5])
def test_powerlaw_matches_analytical(n):
    r = solve_channel(160, 1.0, lambda g: power_law_viscosity(g, K=1.0, n=n))
    ex = powerlaw_velocity(r.y, G=1.0, K=1.0, n=n)
    assert r.converged
    assert l2_error(r.u, ex) < 5e-4


def test_second_order_convergence():
    n = 0.6
    errs = []
    for N in (20, 40, 80, 160):
        r = solve_channel(N, 1.0, lambda g: power_law_viscosity(g, n=n))
        errs.append(l2_error(r.u, powerlaw_velocity(r.y, n=n)))
    p = [observed_order(errs[i], errs[i + 1]) for i in range(3)]
    assert p[-1] > 1.8


def test_wall_shear_stress_balance():
    """Integrated momentum balance requires tau_w = G*h."""
    n = 0.5
    r = solve_channel(200, 1.5, lambda g: power_law_viscosity(g, n=n), h=1.0)
    tw = r.wall_shear_stress(lambda g: power_law_viscosity(g, n=n))
    assert tw == pytest.approx(1.5, rel=2e-3)


def test_bulk_velocity_inversion():
    for n in (1.0, 0.6, 0.3):
        G = pressure_gradient_for_bulk(1.0, K=2.0, n=n)
        assert bulk_velocity_powerlaw(G, K=2.0, n=n) == pytest.approx(1.0)


def test_newtonian_limit_of_powerlaw():
    y = np.linspace(0, 1, 11)
    assert np.allclose(powerlaw_velocity(y, n=1.0), poiseuille_velocity(y))
