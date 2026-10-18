# Numerical Study of Newtonian and Shear-Thinning Channel Flow

**Python / MATLAB · Numerical Methods · Non-Newtonian Fluid Mechanics**

A finite-volume study of steady, pressure-driven flow between parallel plates for a
Newtonian fluid and a shear-thinning fluid described by the **power-law (Ostwald–de Waele)
constitutive relation**, validated against analytical solutions and verified with a
systematic grid-convergence study.

| What the project does | Where |
|---|---|
| Numerical model of pressure-driven channel flow (Newtonian + power-law) | `src/channelflow/solver.py` |
| Velocity profiles, wall shear stress, flow-behaviour changes across rheological parameters (n, K) | `results/figures/fig1, fig2, fig5, fig6`, `docs/report.md` §5 |
| Newtonian case vs analytical Poiseuille flow | `fig3`, `tests/test_solver.py` |
| Numerical consistency: grid refinement, observed order, GCI, momentum balance | `fig4`, `results/table3_*.csv`, `table4_gci.csv` |
| Extension: bounded Carreau fluid (no singular centreline viscosity) | `fig7` |

## Governing problem

Fully developed flow, constant pressure gradient `G = -dp/dx`, half-channel `0 ≤ y ≤ h` by symmetry:

```
d/dy [ mu(|du/dy|) du/dy ] = -G ,   du/dy(0)=0 ,   u(h)=0 ,   mu = K |du/dy|^(n-1)
```

`n = 1` is Newtonian, `n < 1` shear-thinning. Exact solution:
`u(y) = n/(n+1) (G/K)^(1/n) [h^((n+1)/n) − y^((n+1)/n)]`.

## Method

Second-order conservative finite-volume scheme, face viscosity from the face shear rate,
non-linearity solved by under-relaxed Picard iteration, tridiagonal systems solved directly.
See `docs/report.md` for the full derivation.

## Headline results (reproducible with one command)

* Newtonian case reproduces Poiseuille flow to round-off (relative L2 error ≈ 1e-15).
* Power-law solutions converge at the **expected second order** (observed order 2.00–2.05 for
  n = 0.8, 0.6, 0.4); GCI on bulk velocity and wall shear rate < 0.002 % on the 320-cell grid.
* At equal flow rate, wall shear rate rises with shear-thinning (3.00 → 3.25 → 3.67 → 4.50 → 7.00
  for n = 1 → 0.2) and the profile flattens (u_max/U_b: 1.50 → 1.17).
* Integrated momentum balance `τ_w = G·h` is satisfied to ≤ 4e-5 relative error.

## Run it

```bash
pip install -r requirements.txt
python scripts/run_all.py        # regenerates every figure + CSV table (~20 s)
python -m pytest -q              # 8 verification tests
```
MATLAB: open `matlab/channel_flow.m` and run (port of the Python solver).

## Repository layout

```
src/channelflow/   rheology.py  analytical.py  solver.py  metrics.py
scripts/run_all.py reproduces all results
tests/             unit + verification tests
matlab/            MATLAB port
results/           figures/ and CSV tables
docs/report.md     full write-up (theory, method, results, limitations)
```

## Limitations

1-D fully developed flow only (no entrance region, no inertia); inelastic generalised-Newtonian
rheology (no viscoelasticity or yield stress); the pure power law has an unphysical infinite
zero-shear viscosity, handled by a tiny regularisation (and removed entirely by the Carreau extension).
See `docs/report.md` §7.

## References

1. R. B. Bird, R. C. Armstrong, O. Hassager, *Dynamics of Polymeric Liquids, Vol. 1* (Wiley).
2. R. P. Chhabra, J. F. Richardson, *Non-Newtonian Flow and Applied Rheology* (Butterworth-Heinemann).
3. P. J. Roache, "Perspective: A method for uniform reporting of grid refinement studies," *J. Fluids Eng.* 116 (1994) 405–413.
4. P. J. Carreau, "Rheological equations from molecular network theories," *Trans. Soc. Rheol.* 16 (1972) 99–127.

## License
MIT
