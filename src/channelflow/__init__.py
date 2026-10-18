"""Numerical study of Newtonian and shear-thinning (power-law) channel flow."""
from .rheology import power_law_viscosity, carreau_viscosity, newtonian_viscosity
from .analytical import (
    poiseuille_velocity, powerlaw_velocity, bulk_velocity_powerlaw,
    pressure_gradient_for_bulk, wall_shear_rate_powerlaw,
)
from .solver import solve_channel, ChannelResult
from .metrics import l2_error, linf_error, observed_order, gci

__all__ = [
    "power_law_viscosity", "carreau_viscosity", "newtonian_viscosity",
    "poiseuille_velocity", "powerlaw_velocity", "bulk_velocity_powerlaw",
    "pressure_gradient_for_bulk", "wall_shear_rate_powerlaw",
    "solve_channel", "ChannelResult",
    "l2_error", "linf_error", "observed_order", "gci",
]
