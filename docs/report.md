# Numerical Study of Newtonian and Shear-Thinning Channel Flow

## 1. Objective
Develop a numerical model of pressure-driven channel flow that compares Newtonian and
shear-thinning behaviour using a power-law constitutive relation; examine velocity profiles,
wall shear stress and flow-behaviour changes across rheological parameters; and compare the
Newtonian case with analytical Poiseuille flow while assessing numerical consistency through
grid-convergence checks.

## 2. Physical model
Steady, laminar, incompressible, fully developed flow between parallel plates at `y = ±h`, driven
by a constant pressure gradient `G = −dp/dx`. The x-momentum equation reduces to

    dτ/dy = −G   ⇒   τ(y) = G·y   (symmetry at y = 0),  τ_w = G·h

This stress distribution is rheology-independent, so rheology only decides how velocity responds.

**Constitutive models** (generalised Newtonian, τ = μ_app(γ̇)·γ̇):

* Newtonian: μ = const
* Power law: μ_app = K γ̇^(n−1); n = 1 Newtonian, n < 1 shear-thinning
* Carreau (extension): μ_app = μ∞ + (μ0 − μ∞)[1 + (λγ̇)²]^((n−1)/2)

## 3. Analytical benchmarks
* Newtonian (Poiseuille): `u = G/(2μ)(h² − y²)`, `U_b = G h²/(3μ)`, `γ̇_w = 3U_b/h`.
* Power law: `u = n/(n+1)(G/K)^(1/n)[h^((n+1)/n) − |y|^((n+1)/n)]`,
  `U_b = n/(2n+1)(G/K)^(1/n) h^((n+1)/n)`, `γ̇_w = (2n+1)/n · U_b/h`, `u_max/U_b = (2n+1)/(n+1)`.

## 4. Numerical method
Half channel, `N` uniform intervals, unknowns `u_0 … u_{N−1}`, `u_N = 0`. Conservative scheme

    [ μ_{i+½}(u_{i+1}−u_i) − μ_{i−½}(u_i−u_{i−1}) ] / Δy² = −G ,   μ_{i+½} = μ(|u_{i+1}−u_i|/Δy)

with a half-cell symmetry condition at `y = 0` (`2μ_{½}(u_1−u_0)/Δy² = −G`). The non-linearity is solved
by Picard iteration (μ lagged, under-relaxation 0.5, stop at relative update < 1e-10), each linear
system is tridiagonal and solved by banded LU. The power-law viscosity diverges at the centreline for n<1,
so a regularisation `γ̇² → γ̇² + ε²` with ε = 1e-8 is used; its effect is far below discretisation error.
Wall shear rate uses a second-order one-sided stencil; bulk velocity uses the trapezoidal rule.

Comparison fairness: shear-thinning is compared at **equal bulk velocity** (G adjusted analytically per n)
so profiles are comparable, and also at **equal pressure gradient** (Table 2).

## 5. Results (all numbers from `scripts/run_all.py`)

### 5.1 Newtonian validation (Fig. 3)
Relative L2 error vs Poiseuille ≈ 7e-16 at N = 40. This is expected: a second-order scheme differentiates
a parabola exactly, so agreement to round-off confirms correct assembly, boundary conditions and symmetry treatment.

### 5.2 Velocity profiles and wall shear (Fig. 1, Table 1; U_b = 1, K = h = 1, N = 200)

| n | G needed | u_max/U_b | wall shear rate γ̇_w | wall shear stress τ_w | Picard its |
|---|---|---|---|---|---|
| 1.0 | 3.000 | 1.500 | 3.000 | 3.000 | 1 |
| 0.8 | 2.567 | 1.444 | 3.250 | 2.567 | 41 |
| 0.6 | 2.181 | 1.375 | 3.667 | 2.181 | 59 |
| 0.4 | 1.825 | 1.286 | 4.500 | 1.825 | 94 |
| 0.2 | 1.476 | 1.167 | 7.000 | 1.476 | 421 |

* Shear-thinning flattens the profile (plug-like core): at the centre the shear rate is ≈ 0 and viscosity is
  highest, while the shear is concentrated near the walls where viscosity is lowest (Fig. 1, right panels).
* Wall shear **rate** rises as n falls (3 → 7), but wall shear **stress** *falls* (3.00 → 1.48) because the
  pressure gradient needed for the same flow decreases — thinning near the wall lowers resistance.
* Numerical vs analytical: L2 error ≤ 3e-5 for all n; wall-stress balance `τ_w = G h` holds to ≤ 3.7e-5.

### 5.3 Equal pressure gradient (Fig. 2, Table 2; G = K = h = 1)
Because `μ_app(γ̇ = 1) = K` for every n and τ_w = 1 here, all fluids share the same viscosity at the wall
shear rate; in the interior the shear rate is lower, so viscosity is *higher* and `U_b/U_b,Newt` =
0.923, 0.818, 0.667, 0.429 for n = 0.8 … 0.2. The sign of the flow-rate difference therefore depends on the
reference chosen: since `U_b ∝ (G/K)^(1/n)`, shear-thinning fluids carry *more* flow than the matched
Newtonian fluid once `G > K` (high driving pressure), and less when `G < K`. This is a point to state
explicitly when comparing fluids with different n.

### 5.4 Effect of K (Fig. 5 right)
At fixed G, `U_b ∝ K^(−1/n)`; the numerical points sit on the analytical line. Fig. 6 maps τ_w over (n, K).

### 5.5 Grid convergence (Fig. 4, Tables 3–4)

| n | order p (N = 320→640) | L2 error N=40 | L2 error N=640 |
|---|---|---|---|
| 0.8 | 2.005 | 2.6e-5 | 9.9e-8 |
| 0.6 | 2.002 | 7.2e-5 | 2.8e-7 |
| 0.4 | 1.999 | 1.9e-4 | 7.2e-7 |

* The observed order is 2.00–2.05 on every refinement, matching the formal second-order accuracy, for velocity
  and wall shear rate.
* Roache GCI (N = 80/160/320): apparent order ≈ 2.0, and GCI of the fine grid for U_b and γ̇_w is
  ≤ 0.002 %; Richardson-extrapolated values match analytical ones.
* The n = 1 case is excluded from these plots because the exact parabola leaves only round-off error.

### 5.6 Solver behaviour
Picard iterations grow as n decreases (41 → 421 for n = 0.8 → 0.2 at N = 200, and ~360 at N = 640 for n = 0.4):
stronger non-linearity and the near-singular centreline viscosity slow successive substitution. A Newton
linearisation would cut this and is an obvious improvement.

### 5.7 Carreau extension (Fig. 7)
With a bounded viscosity (μ0 = 1, μ∞ = 0.01, n = 0.4), small λ gives Newtonian-plateau (parabolic) flow and large
λ approaches power-law plug-like behaviour, with the singular centreline viscosity replaced by the finite μ0.
Only numerical results are shown (no closed form).

## 6. Conclusions
1. The finite-volume/Picard model reproduces Poiseuille flow exactly and the power-law analytical solution at
   second-order accuracy.
2. Shear-thinning flattens the velocity profile and concentrates shear at the walls; at equal flow rate the wall
   shear rate increases while the required pressure gradient and wall shear stress decrease.
3. The consistency checks (analytical comparison, order of accuracy, GCI, `τ_w = G h` balance) all agree.

## 7. Limitations and future work
* 1-D, fully developed, no inertia; no entrance length or 2-D/3-D effects.
* Purely viscous (inelastic) model; no viscoelasticity, thixotropy or yield stress (Herschel–Bulkley, Bingham are natural next steps).
* Regularisation at the centreline for the pure power law (ε = 1e-8); Carreau avoids this.
* Picard iteration slows for strong shear-thinning; Newton or Anderson acceleration would help.
* Extensions: pipe (axisymmetric) geometry, non-uniform (wall-clustered) grids, 2-D channel with a CFD solver, fitting n and K to real polymer rheometry data.

## 8. References
See `README.md`.
