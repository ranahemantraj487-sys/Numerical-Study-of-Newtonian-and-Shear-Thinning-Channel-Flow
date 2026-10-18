"""Reproduce every figure and table in the report.   Usage:  python scripts/run_all.py"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from channelflow import *

FIG = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
TAB = ROOT / "results"
plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

H, K = 1.0, 1.0
N_FINE = 200
NS = [1.0, 0.8, 0.6, 0.4, 0.2]
col = dict(zip(NS, plt.cm.viridis(np.linspace(0, 0.9, len(NS)))))
pl = lambda n, K=K: (lambda g: power_law_viscosity(g, K=K, n=n))

# ---------------------------------------------------------------- 1. profiles at fixed bulk velocity
rows = []
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
for n in NS:
    G = pressure_gradient_for_bulk(1.0, K=K, n=n, h=H)
    r = solve_channel(N_FINE, G, pl(n), h=H)
    ex = powerlaw_velocity(r.y, G=G, K=K, n=n, h=H)
    ax[0].plot(r.u / r.bulk_velocity, r.y, color=col[n], lw=2, label=f"n = {n}")
    ax[0].plot(ex / 1.0, r.y, "k--", lw=0.8, alpha=0.7)
    gd = r.shear_rate()
    ax[1].plot(r.y, gd, color=col[n], lw=2, label=f"n = {n}")
    ax[2].plot(r.y, pl(n)(gd), color=col[n], lw=2, label=f"n = {n}")
    tw = r.wall_shear_stress(pl(n))
    rows.append(dict(n=n, G_required=G, Ub_num=r.bulk_velocity,
                     umax_over_Ub=r.u[0] / r.bulk_velocity,
                     gamma_w_num=r.wall_shear_rate(),
                     gamma_w_exact=wall_shear_rate_powerlaw(1.0, n, H),
                     tau_w_num=tw, tau_w_exact=G * H,
                     tau_w_rel_err=abs(tw - G * H) / (G * H),
                     L2_error=l2_error(r.u, ex), Linf_error=linf_error(r.u, ex),
                     picard_iters=r.iterations))
ax[0].set(xlabel=r"$u/U_b$", ylabel="y/h", title="Velocity profiles at equal flow rate (dashed: analytical)")
ax[1].set(xlabel="y/h", ylabel=r"$|du/dy|$", title="Shear-rate distribution")
ax[2].set(xlabel="y/h", ylabel=r"$\mu_{app}$", title="Apparent viscosity", yscale="log")
for a in ax: a.legend()
fig.tight_layout(); fig.savefig(FIG / "fig1_profiles_fixed_flowrate.png"); plt.close(fig)
t1 = pd.DataFrame(rows); t1.to_csv(TAB / "table1_fixed_flowrate.csv", index=False)

# ---------------------------------------------------------------- 2. equal pressure gradient
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
rows = []
for n in NS:
    r = solve_channel(N_FINE, 1.0, pl(n), h=H)
    ex = powerlaw_velocity(r.y, G=1.0, K=K, n=n, h=H)
    ax[0].plot(r.u, r.y, color=col[n], lw=2, label=f"n = {n}")
    rows.append(dict(n=n, umax=r.u[0], Ub=r.bulk_velocity, Ub_exact=bulk_velocity_powerlaw(1.0, K, n, H),
                     Ub_over_Newtonian=r.bulk_velocity / (1 / 3.0), L2_error=l2_error(r.u, ex)))
t2 = pd.DataFrame(rows); t2.to_csv(TAB / "table2_fixed_pressure_gradient.csv", index=False)
ns = np.linspace(0.2, 1.0, 60)
ax[1].plot(ns, [bulk_velocity_powerlaw(1.0, K, n, H) * 3 for n in ns], "k-", label="analytical")
ax[1].plot(t2.n, t2.Ub_over_Newtonian, "ro", label="numerical")
ax[0].set(xlabel="u", ylabel="y/h", title=r"Equal $-dp/dx$ ($G=K=h=1$)"); ax[0].legend()
ax[1].set(xlabel="power-law index n", ylabel=r"$U_b / U_{b,Newtonian}$", title=r"Bulk velocity vs Newtonian ($\mu=K$), $G=K=h=1$"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "fig2_equal_pressure_gradient.png"); plt.close(fig)

# ---------------------------------------------------------------- 3. Newtonian vs Poiseuille
r = solve_channel(40, 1.0, lambda g: newtonian_viscosity(g, 1.0), h=H)
ex = poiseuille_velocity(r.y, 1.0, 1.0, H)
newt_err = l2_error(r.u, ex)
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ex, r.y, "k-", lw=2, label="analytical Poiseuille"); ax[0].plot(r.u, r.y, "ro", ms=4, label="numerical (N=40)")
ax[0].set(xlabel="u", ylabel="y/h", title="Newtonian validation"); ax[0].legend()
ax[1].plot(r.y, np.abs(r.u - ex), "b-"); ax[1].set(xlabel="y/h", ylabel="|error|", title=f"Pointwise error (rel. L2 = {newt_err:.1e})")
fig.tight_layout(); fig.savefig(FIG / "fig3_newtonian_vs_poiseuille.png"); plt.close(fig)

# ---------------------------------------------------------------- 4. grid convergence
Ns = [10, 20, 40, 80, 160, 320, 640]
rows = []
for n in [1.0, 0.8, 0.6, 0.4]:
    G = 1.0
    ex_gw = (G / K) ** (1 / n) * 1.0  # |du/dy|_w = (tau_w/K)^(1/n)
    prev = None
    for N in Ns:
        r = solve_channel(N, G, pl(n), h=H)
        ex = powerlaw_velocity(r.y, G=G, K=K, n=n, h=H)
        e = l2_error(r.u, ex)
        rows.append(dict(n=n, N=N, dy=H / N, L2_error=e, Linf_error=linf_error(r.u, ex),
                         order=observed_order(prev, e) if prev else np.nan,
                         Ub=r.bulk_velocity, Ub_exact=bulk_velocity_powerlaw(G, K, n, H),
                         gamma_w=r.wall_shear_rate(), gamma_w_exact=ex_gw,
                         gamma_w_rel_err=abs(r.wall_shear_rate() - ex_gw) / ex_gw,
                         iters=r.iterations))
        prev = e
gc = pd.DataFrame(rows); gc.to_csv(TAB / "table3_grid_convergence.csv", index=False)

# GCI on bulk velocity and wall shear rate using the three finest grids (160, 320, 640 -> use 80,160,320 too)
g_rows = []
for n in [0.8, 0.6, 0.4]:   # n=1 omitted: parabola is reproduced to round-off, so GCI is undefined
    s = gc[gc.n == n].set_index("N")
    for q in ["Ub", "gamma_w"]:
        p, fext, g = gci(s.loc[320, q], s.loc[160, q], s.loc[80, q])
        g_rows.append(dict(n=n, quantity=q, f_N80=s.loc[80, q], f_N160=s.loc[160, q], f_N320=s.loc[320, q],
                           p_apparent=p, f_extrapolated=fext, GCI_fine_percent=100 * g))
gci_df = pd.DataFrame(g_rows); gci_df.to_csv(TAB / "table4_gci.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for n in [0.8, 0.6, 0.4]:   # Newtonian case is exact to round-off (see Fig. 3)
    s = gc[gc.n == n]
    ax[0].loglog(s.dy, s.L2_error.clip(lower=1e-17), "o-", label=f"n = {n}")
    ax[1].loglog(s.dy, s.gamma_w_rel_err.clip(lower=1e-17), "s-", label=f"n = {n}")
d = np.array([H / 640, H / 10])
for a in ax: a.loglog(d, 0.5 * d**2, "k--", label=r"slope 2")
ax[0].set(xlabel=r"$\Delta y/h$", ylabel="relative L2 error (velocity)", title="Grid convergence: velocity")
ax[1].set(xlabel=r"$\Delta y/h$", ylabel="relative error", title="Grid convergence: wall shear rate")
for a in ax: a.legend()
fig.tight_layout(); fig.savefig(FIG / "fig4_grid_convergence.png"); plt.close(fig)

# ---------------------------------------------------------------- 5. consistency: Picard convergence + K sweep
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for n in [0.8, 0.6, 0.4, 0.2]:
    r = solve_channel(100, 1.0, pl(n), h=H)
    ax[0].semilogy(r.residual_history, label=f"n = {n}  ({r.iterations} its)")
ax[0].set(xlabel="Picard iteration", ylabel="relative update", title="Non-linear solver convergence"); ax[0].legend()
Ks = np.array([0.25, 0.5, 1, 2, 4]); n_k = 0.6
Ub_k = [solve_channel(100, 1.0, pl(n_k, K=k), h=H).bulk_velocity for k in Ks]
ax[1].loglog(Ks, Ub_k, "ro", label="numerical")
ax[1].loglog(Ks, [bulk_velocity_powerlaw(1.0, k, n_k, H) for k in Ks], "k-", label=r"analytical $\propto K^{-1/n}$")
ax[1].set(xlabel="consistency index K", ylabel=r"$U_b$", title=f"Effect of K at fixed $-dp/dx$ (n={n_k})"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "fig5_solver_and_K_sweep.png"); plt.close(fig)

# ---------------------------------------------------------------- 6. n-K parameter map of wall shear stress at fixed flow
fig, ax = plt.subplots(figsize=(5.5, 4.4))
nn = np.linspace(0.3, 1.0, 40); KK = np.logspace(-1, 1, 40)
TW = np.array([[KK_ * wall_shear_rate_powerlaw(1.0, n_, H) ** n_ for n_ in nn] for KK_ in KK])
cs = ax.contourf(nn, KK, np.log10(TW), 20, cmap="magma"); ax.set_yscale("log")
plt.colorbar(cs, label=r"$\log_{10}\,\tau_w$"); ax.set(xlabel="n", ylabel="K", title=r"Wall shear stress at $U_b=h=1$")
fig.tight_layout(); fig.savefig(FIG / "fig6_wall_shear_map.png"); plt.close(fig)

# ---------------------------------------------------------------- 7. Carreau extension (numerical only)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
G = 1.0
carr = lambda lam: (lambda g: carreau_viscosity(g, mu0=1.0, mu_inf=0.01, lam=lam, n=0.4))
for lam, c in zip([0.01, 1, 10, 100], plt.cm.plasma(np.linspace(0.1, 0.85, 4))):
    r = solve_channel(N_FINE, G, carr(lam), h=H, relax=0.5)
    ax[0].plot(r.u, r.y, color=c, lw=2, label=fr"$\lambda$ = {lam}")
    ax[1].plot(r.y, carr(lam)(r.shear_rate()), color=c, lw=2, label=fr"$\lambda$ = {lam}")
ax[0].set(xlabel="u", ylabel="y/h", title="Carreau fluid ($\\mu_0=1,\\ \\mu_\\infty=0.01,\\ n=0.4$)"); ax[0].legend()
ax[1].set(xlabel="y/h", ylabel=r"$\mu_{app}$", yscale="log", title="Apparent viscosity"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "fig7_carreau_extension.png"); plt.close(fig)

# ---------------------------------------------------------------- print summary
pd.set_option("display.width", 200, "display.float_format", lambda x: f"{x:.4g}")
print("\nTABLE 1 fixed flow rate\n", t1.to_string(index=False))
print("\nTABLE 2 fixed pressure gradient\n", t2.to_string(index=False))
print(f"\nNewtonian (N=40) rel L2 error vs Poiseuille: {newt_err:.3e}")
print("\nTABLE 3 grid convergence\n", gc.to_string(index=False))
print("\nTABLE 4 GCI\n", gci_df.to_string(index=False))
