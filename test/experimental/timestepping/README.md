# Benchmark Time-Stepping Methods {#benchmark-timestepping-methods}

## Verification of Multi-Stage Time Integration Schemes

**Problem Description**

DuMux implements a family of one-step multi-stage time integration schemes
in Shu-Osher form @cite ShuOsher1988 (see `MultiStageMethod`). Given the
semi-discrete problem

$$
\frac{\partial M(x,t)}{\partial t} - R(x,t) = 0,
$$

an $m$-stage scheme advances the solution from $t^n$ to $t^{n+1}=t^n+\Delta t^n$ via

$$
\sum_{k=0}^{i}\Big[\alpha_{ik}\,M\big(x^{(k)}, t^n+d_k\Delta t^n\big)
+ \beta_{ik}\,\Delta t^n\, R\big(x^{(k)}, t^n+d_k\Delta t^n\big)\Big] = 0,
\qquad i=1,\ldots,m,
$$

with $x^{(0)}=x^n$ and $x^{n+1}=x^{(m)}$. Restricting the inner sum to $k\le i$
yields explicit schemes (if $\beta_{ii}=0$) or diagonally implicit (DIRK)
schemes (if $\beta_{ii}>0$). The following schemes are implemented and
verified here:

| Scheme | id | Stages $m$ | Order $p$ | Type |
|---|---|---|---|---|
| Explicit Euler | `explicit_euler` | 1 | 1 | explicit |
| Implicit Euler | `implicit_euler` | 1 | 1 | implicit, L-stable |
| Crank-Nicolson | `crank_nicolson` | 1 | 2 | implicit, A-stable |
| Heun (explicit RK2) | `heun` | 2 | 2 | explicit |
| DIRK2 (Alexander) | `dirk2` | 2 | 2 | implicit, L-stable |
| Explicit RK3 | `rk3` | 3 | 3 | explicit |
| DIRK3 (Alexander) | `dirk3` | 3 | 3 | implicit, L-stable |
| Explicit RK4 | `rk4` | 4 | 4 | explicit |

Crank-Nicolson is the theta scheme with $\theta=0.5$, and the two DIRK
schemes are taken from Alexander (1977) @cite Alexander1977.

Two complementary, purely time-discrete verification tests are used (the
spatial operator is trivial in both, so only the temporal error is
measured):

1. **Observed order of convergence**: a scalar initial value problem is
   integrated to $t=1$ with successively halved step sizes $\Delta t$, and
   the error is compared to its expected asymptotic decay rate
   $\mathcal{O}(\Delta t^p)$ (two different right-hand sides are used, see
   *Analytical Solution* below).
2. **Linear stability region**: for the Dahlquist test equation
   $\dot u = \lambda u$, one step gives $u^{n+1}=R(z)\,u^n$ with
   $z=\lambda\Delta t \in \mathbb{C}$ and stability function $R$. The
   region $\{z : |R(z)|\le 1\}$ is traced numerically and compared to the
   analytical $R(z)$.

**Analytical Solution**

Two scalar test equations are used.

For the single-step accuracy checks, $\dot u = e^t,\ u(0)=0$ has the exact
solution $u(t) = e^t-1$.

For the convergence-sweep figure below, $\dot u = -u + e^t,\ u(0)=0$ has the
exact solution $u(t) = \sinh(t)$. Unlike the first equation, its right-hand
side depends on $u$ (see the note in *Results* on why this matters).

For an explicit Runge-Kutta scheme of classical order $p$ with $p$ stages,
the stability function is the truncated exponential

$$
R(z) = \sum_{k=0}^{p}\frac{z^k}{k!},
$$

i.e. the order-$p$ Taylor polynomial of $e^z$. Its boundary $|R(z)|=1$ is a
closed curve; the *negative real-axis cutoff* is the largest $|\mathrm{Re}(z)|$
for which $z<0$ lies on this boundary (the maximum stable step size for a
decaying mode with eigenvalue $\lambda=z/\Delta t<0$).

For implicit schemes, **A-stability** requires $|R(z)|\le 1$ for all
$\mathrm{Re}(z)\le 0$ (unconditional stability for decaying modes), and
**L-stability** additionally requires $R(z)\to 0$ as $\mathrm{Re}(z)\to-\infty$
(strong damping of stiff modes).

**Setup**

Both tests advance a single scalar degree of freedom with the
`MultiStageTimeStepper` and require neither a mesh nor a linear solver
beyond a $1\times1$ system:

- `test_timestepmethods.cc` solves $\dot u = e^t$ with all 8 schemes,
  checks the Shu-Osher coefficients for consistency and the single-step
  error at $\Delta t=0.01$ against per-scheme tolerances. For the
  convergence-sweep figure, it additionally solves $\dot u = -u + e^t$,
  $u(0)=0$ with all 8 schemes for
  $\Delta t \in \{0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625\}$ and writes the
  errors at $t=1$ to `test_timestepmethods_data.json`.
- `test_timestepmethods_stabilityregions.cc` evaluates the stability
  function $R(z)$ of each scheme directly from its Shu-Osher coefficients,
  traces the boundary $|R(z)|=1$ for the explicit schemes by continuation
  (parametrizing $R(z(\theta))=e^{i\theta}$ and following the curve with
  Newton's method as $\theta$ increases), checks A-/L-stability of the
  implicit schemes on a grid in the left half-plane, and writes
  `test_timestepmethods_stabilityregions_data.json`.


**Benchmark: observed order of convergence**

![Observed order of temporal convergence](timestepping_convergence.png)

Each curve shows the absolute error $|u_{\Delta t}(1) - u_{\text{exact}}(1)|$
versus $\Delta t$ on a log-log scale for $\dot u = -u + e^t,\ u(0)=0$
(exact solution $u(t)=\sinh(t)$); the triangles annotate reference slopes
$\mathcal{O}(\Delta t^p)$ for $p=1,\ldots,4$ (rise $p$ over run $1$). All
schemes attain their expected asymptotic order, or better.

Note: for the simpler equation $\dot u = e^t$ used in the single-step
checks above, the right-hand side is independent of $u$, so each scheme
reduces to a quadrature rule for $\int e^t\,dt$: Crank-Nicolson and Heun
both reduce to the trapezoidal rule and are therefore identical, and the
nodes and weights of the explicit RK3 scheme coincide with Simpson's rule
(as do RK4's), so both converge with fourth order for that equation. The
convergence-sweep equation $\dot u=-u+e^t$ above avoids this degeneracy,
since its right-hand side depends on $u$.

**Benchmark: Linear stability regions**

| Scheme | negative real-axis cutoff |
|---|---|
| Explicit Euler | 2.0 |
| Heun (RK2) | 2.0 |
| Explicit RK3 | 2.5127 |
| Explicit RK4 | 2.7853 |

![Stability regions of explicit schemes](timestepping_stability_explicit.png)

![Numerical vs. analytical stability boundaries](timestepping_stability_comparison.png)

The comparison plot overlays the numerically traced boundary (solid line)
with the analytical truncated-exponential boundary $R(z)=\sum_{k=0}^p z^k/k!$
(crosses); the two coincide to within the Newton tolerance.

![Stability regions of implicit schemes](timestepping_stability_implicit.png)

For the implicit schemes, implicit Euler, DIRK2, and DIRK3 are L-stable (the
stable region extends unboundedly into the left half-plane, with
$R(z)\to0$), while Crank-Nicolson is A-stable but not L-stable
($|R(z)|\to1$ as $\mathrm{Re}(z)\to-\infty$).

**Reproducing the figures**

After building the `test_timestepmethods` and
`test_timestepmethods_stabilityregions` executables, run:

```bash
python3 run_benchmarks.py <path-to-build-dir>
```

This runs both executables to (re-)generate the JSON data and produces all
four figures above via `plot_timestepmethods_convergence.py` and
`plot_timestepmethods_stabilityregions.py`. If run from within
`<build-dir>/test/experimental/timestepping`, the build directory argument
can be omitted. Pass `--copy-to-doc` to additionally copy the figures into
`doc/doxygen/images/` (used to update the images shown on this page).
