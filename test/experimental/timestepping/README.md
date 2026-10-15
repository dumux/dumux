# Benchmark Time-Stepping Methods {#benchmark-timestepping-methods}

## Verification of Multi-Stage Time Integration Schemes

**Problem Description**

DuMux implements a family of one-step multi-stage time integration schemes
in Shu-Osher form @cite ShuOsher1988 (see `MultiStageMethod`). Given the
semi-discrete problem

$$
\frac{\partial M(x,t)}{\partial t} + R(x,t) = 0,
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
| Qin-Zhang | `qin_zhang` | 2 | 2 | implicit, A-stable, symplectic |
| Explicit RK3 | `rk3` | 3 | 3 | explicit |
| DIRK3 (Alexander) | `dirk3` | 3 | 3 | implicit, L-stable |
| Explicit RK4 | `rk4` | 4 | 4 | explicit |
| DIRK4 (Hairer-Wanner) | `dirk4` | 5 | 4 | implicit, L-stable |

Crank-Nicolson is the theta scheme with $\theta=0.5$, the DIRK2 and DIRK3
schemes are taken from Alexander (1977) @cite Alexander1977, DIRK4 is the
five-stage stiffly accurate scheme with diagonal coefficient $1/4$ of Hairer
and Wanner (1996) @cite HairerWanner1996 (Section IV.6), and the Qin-Zhang scheme
is the two-stage second-order *symplectic* DIRK of Qin and Zhang (1992)
@cite QinZhang1992. Being symplectic, it conserves quadratic invariants
exactly; this property is verified separately (see *Benchmark: symplecticity*
below).

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

The tests require neither a mesh nor a linear solver beyond a small dense
system:

- `test_timestepmethods.cc` advances a single scalar degree of freedom with
  the `MultiStageTimeStepper` and solves $\dot u = e^t$ with all 10 schemes,
  checks the Shu-Osher coefficients for consistency and the single-step
  error at $\Delta t=0.01$ against per-scheme tolerances. For the
  convergence-sweep figure, it additionally solves $\dot u = -u + e^t$,
  $u(0)=0$ with all 10 schemes for
  $\Delta t \in \{0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625\}$ and writes the
  errors at $t=1$ to `test_timestepmethods_data.json`.
- `test_timestepmethods_stabilityregions.cc` evaluates the stability
  function $R(z)$ of each scheme directly from its Shu-Osher coefficients,
  traces the boundary $|R(z)|=1$ for the explicit schemes by continuation
  (parametrizing $R(z(\theta))=e^{i\theta}$ and following the curve with
  Newton's method as $\theta$ increases), checks A-/L-stability of the
  implicit schemes on a grid in the left half-plane, and writes
  `test_timestepmethods_stabilityregions_data.json`.
- `test_timestepmethods_symplectic.cc` verifies the symplecticity of the
  Qin-Zhang scheme, both algebraically (the Cooper condition on its Butcher
  tableau) and numerically on the $2\times2$ harmonic-oscillator system
  (phase-space area preservation and energy conservation), and writes
  `test_timestepmethods_symplectic_data.json`.


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

For the implicit schemes, implicit Euler, DIRK2, DIRK3 and DIRK4 are L-stable (the
stable region extends unboundedly into the left half-plane, with
$R(z)\to0$), while Crank-Nicolson and the symplectic Qin-Zhang scheme are
A-stable but not L-stable ($|R(z)|\to1$ as $\mathrm{Re}(z)\to-\infty$). For
Qin-Zhang the stability function is $R(z)=\big((4+z)/(4-z)\big)^2$, so
$|R(iy)|=1$ on the imaginary axis: it introduces no artificial damping, which
is the flip side of its symplecticity.

**Benchmark: symplecticity**

A Runge-Kutta method with weights $b_i$ and coefficients $a_{ij}$ is
*symplectic* if and only if it satisfies the condition

$$
b_i a_{ij} + b_j a_{ji} - b_i b_j = 0 \qquad \forall\, i,j,
$$

which is exactly the condition under which the method conserves every
quadratic first integral exactly. `test_timestepmethods_symplectic.cc`
verifies this condition algebraically for the Qin-Zhang tableau (extracted
directly from its Shu-Osher coefficients) and then demonstrates its two
consequences on the harmonic oscillator $\dot q = p,\ \dot p = -q$ (a linear
Hamiltonian system with quadratic energy $H=(q^2+p^2)/2$):

1. **Phase-space area preservation**: the determinant of the one-step
   propagation matrix is $1$ to machine precision for Qin-Zhang, whereas it
   is $<1$ for implicit Euler and DIRK2 (they contract phase space).
2. **Energy conservation**: integrated over 200 periods, Qin-Zhang conserves
   $H$ to a relative drift below $10^{-10}$, while implicit Euler (1st order)
   damps it away almost entirely and the *same second-order* DIRK2 also
   dissipates it — showing that it is symplecticity, not the order of
   accuracy, that conserves energy.

![Symplecticity of the Qin-Zhang scheme](timestepping_symplectic.png)

The left panel shows the phase-space trajectory of the symplectic Qin-Zhang
scheme (a closed circle) against dissipative implicit Euler (spiralling
inward); DIRK2 is omitted from this panel because it is so weakly dissipative
that it is visually indistinguishable from the symplectic orbit at this
resolution. The right panel shows the energy drift $|E(t)/E_0-1|$ over time
on a logarithmic axis, where that small dissipation becomes visible:
Qin-Zhang stays at the machine-precision floor ($\sim 10^{-15}$), whereas
implicit Euler loses essentially all the energy and the *same second-order*
DIRK2 drifts by $\sim 10^{-4}$ — roughly eleven orders of magnitude more than
the symplectic scheme, confirming that it is the symplecticity, not the order
of accuracy, that conserves the energy.

## Verification of Multi-Step (BDF) Time Integration Schemes

**Problem Description**

A linear $k$-step method (see `MultiStepMethod` and `MultiStepTimeStepper`)
computes the solution $x^{n+1}$ at $t^{n+1}=t^n+\Delta t^n$ from the solutions
at the $k$ previous time levels by solving

$$
\sum_{j=0}^{k}\Big[\alpha_j\,M\big(x^{n+1-j}, t^{n+1-j}\big)
+ \beta_j\,\Delta t^n\, R\big(x^{n+1-j}, t^{n+1-j}\big)\Big] = 0 .
$$

The backward differentiation formula (BDF) of order $k$ evaluates the spatial
operator at the new time level only, $\beta=(1,0,\ldots,0)$, and
differentiates the storage term with the interpolation polynomial through the
time levels $t^{n+1},\ldots,t^{n+1-k}$:
$\alpha_j = \Delta t^n\,\ell_j'(t^{n+1})$, where $\ell_j$ are the Lagrange
basis polynomials of the time levels. The weights therefore depend on the step
sizes $\Delta t^n,\ldots,\Delta t^{n+1-k}$. Because $R$ enters at the new time
level only, equations without a storage term (constraints) hold after every
step.

| Scheme | Order $p$ | Linear stability (constant $\Delta t$) | Used for |
|---|---|---|---|
| BDF1 (implicit Euler) | 1 | A-stable, L-stable | all step sizes |
| BDF2 | 2 | A-stable | step size ratio $\omega=\Delta t^n/\Delta t^{n-1} \le 2$ |
| BDF3 | 3 | $A(86.03^\circ)$ | equal step sizes |
| BDF4 | 4 | $A(73.35^\circ)$ | equal step sizes |
| BDF5 | 5 | $A(51.84^\circ)$ | equal step sizes |

Variable-step BDF2 is zero-stable for $\omega < 1+\sqrt{2}$
@cite Grigorieff1983; the default limit of the implementation is $\omega \le 2$
(adjustable). A step whose step sizes do not satisfy the condition of the
order is taken with the highest lower order whose condition holds (BDF2,
otherwise BDF1). The $A(\alpha)$ angles are those of Hairer and Wanner (1996)
@cite HairerWanner1996 (Section V.2).

A BDF method of order $k$ needs $k$ previous time levels. While fewer are
available, i.e. in the first steps and after `resetHistory()`, the stepper
takes the steps with a startup method: by default DIRK4 for BDF3 to BDF5, and
the lower BDF orders for BDF1 and BDF2. The global error contains the local
error of these steps, so a startup method of order $p-1$ is required for
global order $p$.

**Analytical Solution**

The convergence tests solve $\dot u = -u + e^t$, $u(0)=1$, with the exact
solution $u(t)=\cosh(t)$. Since $\ddot u(0) = 1 \neq 0$, the local error of
the first steps does not vanish, and the order of the startup steps shows in
the global error. The constrained system

$$
\dot u_0 = -u_0 + u_1, \qquad 0 = \sin(t) - u_1, \qquad u_0(0)=0,
$$

has no storage term in the second equation and the exact solution
$u_0(t)=\tfrac12\left(\sin t - \cos t + e^{-t}\right)$, $u_1(t)=\sin t$.

For the Dahlquist test equation $\dot u=\lambda u$ ($M(u)=u$,
$R(u)=-\lambda u$) with constant step size, BDF$k$ gives the recurrence
$\sum_j \alpha_j u^{n+1-j} = z\,u^{n+1}$ with $z=\lambda\Delta t$. It is stable
if all roots $\zeta$ of $\sum_j \alpha_j \zeta^{k-j} - z\,\zeta^k$ lie in the
closed unit disk. The boundary of the stability region is the root locus
$z(\theta)=\sum_j \alpha_j e^{-ij\theta}$, $\theta\in[0,2\pi]$. For variable-step
BDF2 with a constant ratio $\omega$, the equation $\dot u = 0$ gives the roots
$\zeta_1 = 1$ and $\zeta_2 = \omega^2/(1+2\omega)$.

**Setup**

- `test_multistepmethods.cc` checks that the BDF weights reproduce the
  constant-step tables, differentiate polynomials of degree $\le k$ exactly for
  random step sizes and match the closed form of variable-step BDF2; checks the
  number of time levels used for given step sizes; computes the errors at
  $t=1$ for $n=10,20,40,80,160$ steps with constant step sizes, with step sizes
  $t_i = g(i/n)$, $g(s)=s+\tfrac{1}{4\pi}\sin(2\pi s)$ (ratios between
  $0.5$ and $1.5$ of the mean), and with one step shortened to a tenth; checks
  that a step retried with half the step size after a solver failure gives the
  same result as taking the reduced step directly; and solves the constrained
  system. It writes `test_multistepmethods_data.json`.
- `test_multistepmethods_stability.cc` computes the roots of the
  characteristic polynomial (Durand-Kerner iteration), checks that they have
  modulus one on the boundary locus, computes the $A(\alpha)$ angles from the
  locus and checks stability inside the sectors for $|z|$ up to $10^4$, the
  damping of stiff modes and the parasitic root of BDF2. It writes
  `test_multistepmethods_stability_data.json`.

**Benchmark: observed order of convergence**

![Observed order of convergence of the BDF methods](timestepping_bdf_convergence.png)

With constant step sizes, the observed orders between the two finest step sizes
are 0.998, 2.000, 2.983 and 3.972 for BDF1 to BDF4, and 4.94 for BDF5 between
$n=40$ and $n=80$ (for $n=160$ its error approaches the round-off level).
BDF3 started with BDF1 and BDF2 instead of DIRK4 converges with order 2.008:
the local error of the first BDF1 step, which is of order $\Delta t^2$, limits
the global error to second order.

With variable step sizes, BDF1 and BDF2 converge with orders 0.997 and 1.998;
BDF3 takes the steps with BDF2 and converges with order 1.969. The step after
the shortened step has the ratio $\omega=10$ and is taken with BDF1; BDF2 keeps
the order 2.004.

For the constrained system, $|u_1 - \sin t|$ is zero after every step for all
orders, and $u_0$ converges with the orders 0.995, 2.020, 2.973, 3.982 and 4.970
(BDF5 between $n=40$ and $n=80$).

**Benchmark: linear stability**

![Stability regions of the BDF methods](timestepping_bdf_stability.png)

The left panel shows the boundary loci; each method is unstable inside its
curve. BDF1 and BDF2 are A-stable. The dashed lines bound the sectors
$|\arg(-z)| \le \alpha$ of BDF3 to BDF5, whose angles computed from the loci are
$86.03^\circ$, $73.35^\circ$ and $51.84^\circ$. At $z=-10^6$ the largest root
modulus is $1.0\cdot10^{-6}$, $7.1\cdot10^{-4}$, $7.0\cdot10^{-3}$,
$2.3\cdot10^{-2}$ and $4.8\cdot10^{-2}$ for BDF1 to BDF5. The right panel shows
the parasitic root of variable-step BDF2, which reaches the unit circle at
$\omega=1+\sqrt{2}$ and is $0.8$ at the default limit $\omega=2$.

## Verification of the Generalized-α Method for Second Time Derivatives

**Problem Description**

For the semi-discrete problem with inertia operator $I$, storage operator $M$
and spatial operator $R$,

$$
\frac{\partial^2 I(x,t)}{\partial t^2} + \frac{\partial M(x,t)}{\partial t} + R(x,t) = 0,
$$

the generalized-α method (see `GeneralizedAlphaMethod` and
`GeneralizedAlphaTimeStepper`) advances the values of $I$ and $M$ together with
their first and second time derivatives. For $X \in \{I, M\}$, the Newmark
relations

$$
X^{n+1} = X^n + \Delta t\,\dot X^n + \Delta t^2\Big[\big(\tfrac12-\beta\big)\ddot X^n + \beta\,\ddot X^{n+1}\Big],
\qquad
\dot X^{n+1} = \dot X^n + \Delta t\Big[(1-\gamma)\,\ddot X^n + \gamma\,\ddot X^{n+1}\Big]
$$

express the time derivatives at the new time level by $X^{n+1}=X(x^{n+1},t^{n+1})$,
and the method solves

$$
(1-\alpha_m)\,\ddot I^{n+1} + \alpha_m\,\ddot I^{n}
+ (1-\alpha_f)\,\dot M^{n+1} + \alpha_f\,\dot M^{n}
+ (1-\alpha_f)\,R\big(x^{n+1},t^{n+1}\big) + \alpha_f\,R\big(x^{n},t^{n}\big) = 0
$$

for $x^{n+1}$. For $I$ and $M$ linear in $x$, this is the method of Chung and
Hulbert (1993) @cite ChungHulbert1993 for
$\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{r}(\mathbf{u}) = 0$
with the spatial term averaged between the time levels. The parameter sets are
given by the spectral radius $\rho_\infty$ of the amplification matrix in the
limit $\omega\Delta t\to\infty$, where $\omega$ is the angular frequency of an
undamped linear oscillator:

| Scheme | $\alpha_m$ | $\alpha_f$ | Range of $\rho_\infty$ |
|---|---|---|---|
| Chung-Hulbert @cite ChungHulbert1993 | $(2\rho_\infty-1)/(\rho_\infty+1)$ | $\rho_\infty/(\rho_\infty+1)$ | $[0, 1]$ |
| Hilber-Hughes-Taylor (HHT) @cite HilberHughesTaylor1977 | $0$ | $(1-\rho_\infty)/(1+\rho_\infty)$ | $[\tfrac12, 1]$ |
| Wood-Bossak-Zienkiewicz (WBZ) @cite WoodBossakZienkiewicz1980 | $(\rho_\infty-1)/(\rho_\infty+1)$ | $0$ | $[0, 1]$ |
| Newmark @cite Newmark1959 | $0$ | $0$ | $1$ for $(\beta,\gamma)=(\tfrac14,\tfrac12)$ |

The first three use $\gamma=\tfrac12-\alpha_m+\alpha_f$, which gives second
order, and $\beta=\tfrac14(1-\alpha_m+\alpha_f)^2$, which together with
$\alpha_m\le\alpha_f\le\tfrac12$ (satisfied in the given ranges of
$\rho_\infty$) gives unconditional stability for linear problems. The Newmark method takes $\beta$ and $\gamma$
as parameters; the default $(\beta,\gamma)=(\tfrac14,\tfrac12)$ is the average
acceleration method.

An equation without inertia and storage terms, $g(x,t)=0$, is solved in the
form $(1-\alpha_f)\,g^{n+1} + \alpha_f\,g^n = 0$, so its residual is multiplied
by $-\alpha_f/(1-\alpha_f)$ in every step: by $-\rho_\infty$ for
Chung-Hulbert, by $-(1-\rho_\infty)/(2\rho_\infty)$ for HHT, and by zero for
WBZ and Newmark. For Chung-Hulbert with $\rho_\infty=1$, an initial violation
of such an equation persists with alternating sign.

The time derivatives of $I$ and $M$ at the initial time are the Jacobians of
the operators times given time derivatives of the solution
(`initialize(vars, t0, velocity, acceleration)`). `initialize(vars, t0)`
starts problems with inertia at rest (vanishing time derivatives of the
operators) and problems without inertia with $\dot M^0=-R(x^0,t^0)$ and
$\ddot M^0=0$. This rate is consistent if the initial solution satisfies the
equations without storage term.

Without an inertia operator, the stepper integrates first-order problems, and
$M$ follows the Newmark relations. For $\dot u=\lambda u$ ($M(u)=u$,
$R(u)=-\lambda u$), the stepper with the Chung-Hulbert parameters has the
eigenvalues of the first-order generalized-α method of Jansen, Whiting and
Hulbert (2000) @cite JansenWhitingHulbert2000 with the same $\rho_\infty$,
which solves
$(1-\alpha^J_m)\,\dot u^n + \alpha^J_m\,\dot u^{n+1} = \lambda\big((1-\alpha^J_f)\,u^n + \alpha^J_f\,u^{n+1}\big)$
with $u^{n+1}=u^n+\Delta t\big[(1-\gamma^J)\,\dot u^n+\gamma^J\,\dot u^{n+1}\big]$,
$\alpha^J_m=(3-\rho_\infty)/(2(1+\rho_\infty))$, $\alpha^J_f=1/(1+\rho_\infty)$
and $\gamma^J=\tfrac12+\alpha^J_m-\alpha^J_f$, and the additional eigenvalue
$-\alpha_f/(1-\alpha_f)=-\rho_\infty$, by which the residual
$\dot u-\lambda u$ is multiplied in every step. A consistent start does not
excite this eigenvalue.

**Analytical Solution**

The convergence tests solve two decoupled damped oscillators
$\ddot u_i + c_i\,\dot u_i + k_i\,u_i = 0$ with $(c_0,k_0)=(0.2,4)$,
$(c_1,k_1)=(0.5,9)$, $u(0)=(1,0.5)$ and $\dot u(0)=(0,1)$. With
$\omega_i=\sqrt{k_i}$, $\zeta_i=c_i/(2\omega_i)$ and
$\omega_{d,i}=\omega_i\sqrt{1-\zeta_i^2}$, the exact solution is

$$
u_i(t) = e^{-\zeta_i\omega_i t}\left[u_i(0)\cos(\omega_{d,i}t)
+ \frac{\dot u_i(0)+\zeta_i\omega_i u_i(0)}{\omega_{d,i}}\sin(\omega_{d,i}t)\right].
$$

The first-order equation $\dot u=-u+e^t$ has the exact solution
$u(t)=(u(0)-\tfrac12)\,e^{-t}+\tfrac12\,e^t$.

For the undamped oscillator $\ddot u + \omega^2 u = 0$, one step maps
$(u^n,\dot u^n,\ddot u^n)$ to $(u^{n+1},\dot u^{n+1},\ddot u^{n+1})$ with the
amplification matrix $\mathbf{A}(\omega\Delta t)$. Its spectral radius $\rho$
is the largest modulus of its eigenvalues. The principal eigenvalues
$r\,e^{\pm i\bar\theta}$ approximate $e^{\pm i\omega\Delta t}$; they define the
algorithmic damping ratio $\bar\xi=-\ln(r)/\bar\theta$ and the numerical period
$\bar T=2\pi\Delta t/\bar\theta$, compared with the period $T=2\pi/\omega$. For
the Newmark method with $(\beta,\gamma)=(\tfrac14,\tfrac12)$, $r=1$ and
$\bar\theta=2\arctan(\omega\Delta t/2)$.

**Setup**

- `test_generalizedalpha.cc` checks the parameters of the four schemes and
  that the stepper, which advances operator values, reproduces the
  generalized-α method written for displacement, velocity and acceleration in
  60 steps with $\Delta t=0.05$ for a system with two degrees of freedom, a
  non-diagonal mass matrix, a non-symmetric damping matrix and the nonlinear
  stiffness term $0.5\,u_i^3$. It computes the errors of the damped
  oscillators at $t=2$ for $n=25,50,100,200,400$ steps and for $n=200,400$
  steps with the step sizes $t_i=2\,g(i/n)$,
  $g(s)=s+\tfrac{1}{4\pi}\sin(2\pi s)$ (ratios between $0.5$ and $1.5$ of the
  mean), and the amplification
  matrix of the stepper for $10^{-2}\le\omega\Delta t\le10^{6}$ (with
  $\Delta t=1$). For $\dot u=\lambda u$ (storage operator only, $\Delta t=1$),
  it computes the amplification matrix of the stepper on 41 radii
  $10^{-2}\le|\lambda|\le10^{6}$ times 37 angles in the closed left
  half-plane, on the negative real axis and at $\lambda=-10^8$, and compares
  the eigenvalues for Chung-Hulbert with those of Jansen-Whiting-Hulbert and
  $-\rho_\infty$ at seven values of $\lambda$. It computes the order for
  $\dot u=-u+e^t$ with $u(0)=1$ and consistent initial rates, and the errors at
  $t=1$ for $u(0)=2$ and $n=10,20,\ldots,640$ steps with the default start and
  with zero initial rates, and for $n=100,200$ steps with the step sizes
  $t_i=g(i/n)$. It checks that a step retried with half the step size after a
  solver failure gives the same result as taking the reduced step directly. It checks the factor $-\alpha_f/(1-\alpha_f)$ for
  the constraint of the constrained system of the BDF benchmark with
  $u_1(0)-\sin(0)=0.1$ and zero initial rates, checks that the default start
  keeps the constraint for $u(0)=(0.5,0)$ in 20 steps with $\Delta t=0.1$, and
  checks that the Newmark parameters reproduce `NewmarkBeta` in 100 steps of
  $\ddot u=-u$ with $(\beta,\gamma)=(0.3,0.6)$. It writes
  `test_generalizedalpha_data.json`.

**Benchmark: convergence and dissipation**

![Convergence and dissipation of the generalized-α methods](timestepping_generalizedalpha.png)

The largest difference between the stepper and the generalized-α method for
displacement, velocity and acceleration is $4.5\cdot10^{-14}$ over the seven
parameter sets. For the damped
oscillators (top left), the observed orders between $n=200$ and $n=400$ are
2.000 for Newmark and Chung-Hulbert with $\rho_\infty=1$ and $0.8$, 2.004,
2.024, 2.003 and 2.004 for Chung-Hulbert with $\rho_\infty=0.5$ and $0$, HHT
and WBZ with $\rho_\infty=0.8$. With the variable step sizes, the observed
orders are the same to three digits, and the errors at $n=400$ are 1.34 to 1.37
times those with constant step sizes. Chung-Hulbert with $\rho_\infty=1$ and the
Newmark method have the same damping ratio and period (to $10^{-9}$ for
$\omega\Delta t\le 10^4$), and their curves coincide.

The spectral radius (top centre) at $\omega\Delta t=10^6$ is 0.80008, 0.50008
and 0.00010 for Chung-Hulbert with $\rho_\infty=0.8$, $0.5$ and $0$, and
0.80000 for HHT and WBZ with $\rho_\infty=0.8$. The largest spectral radius
over all $\omega\Delta t$ is one for all methods up to $6.3\cdot10^{-6}$ for
Chung-Hulbert with $\rho_\infty=1$, whose three eigenvalues approach $-1$ for
$\omega\Delta t\to\infty$; the roots of the characteristic polynomial are
then only accurate to about the cube root of the round-off. At
$\Delta t/T=0.1$, the damping ratio (bottom left) and the relative period
error (bottom centre) are:

| Scheme | $\bar\xi$ | $(\bar T-T)/T$ |
|---|---|---|
| Newmark, $(\beta,\gamma)=(\tfrac14,\tfrac12)$ | 0 | 0.03207 |
| Chung-Hulbert, $\rho_\infty=1$ | 0 | 0.03207 |
| Chung-Hulbert, $\rho_\infty=0.8$ | 0.00016 | 0.03383 |
| Chung-Hulbert, $\rho_\infty=0.5$ | 0.00378 | 0.04672 |
| Chung-Hulbert, $\rho_\infty=0$ | 0.05499 | 0.12321 |
| HHT, $\rho_\infty=0.8$ | 0.00228 | 0.04017 |
| WBZ, $\rho_\infty=0.8$ | 0.00340 | 0.04192 |

For $\rho_\infty=0.8$, the damping ratio of Chung-Hulbert at $\Delta t/T=0.1$
is $7\,\%$ of that of HHT and $5\,\%$ of that of WBZ.

For $\dot u=\lambda u$ (top right), all parameter sets are A-stable: the
largest spectral radius in the closed left half-plane exceeds one by at most
$3.8\cdot10^{-5}$, for Chung-Hulbert with $\rho_\infty=1$, whose eigenvalues
approach $-1$ for $\lambda\Delta t\to-\infty$. At $\lambda=-10^8$, the
spectral radius is 0.80000, 0.50000 and 0.00007 for Chung-Hulbert with
$\rho_\infty=0.8$, $0.5$ and $0$, 0.80000 for HHT and WBZ, and one for Newmark.
For Chung-Hulbert, the eigenvalues agree with those of Jansen-Whiting-Hulbert
and $-\rho_\infty$ to $3.5\cdot10^{-14}$ ($5.9\cdot10^{-8}$ for
$\rho_\infty=1$). The eigenvalue $-\rho_\infty$ does not depend on $\lambda$,
so the spectral radius of Chung-Hulbert is $\rho_\infty$ wherever the
eigenvalues of Jansen-Whiting-Hulbert are smaller.

For $\dot u=-u+e^t$ with $u(0)=1$ and consistent initial rates, the observed
orders between $n=100$ and $n=200$ are 2.000 for Newmark, HHT, WBZ and
Chung-Hulbert with $\rho_\infty=1$ and $0.8$, and 1.999 and 1.990 for
Chung-Hulbert with $\rho_\infty=0.5$ and $0$. For $u(0)=2$ (bottom right), the
default start gives the order 2.000 between $n=320$ and $n=640$ for all
parameter sets. With zero initial rates, the order is 1.00 for all parameter
sets except Chung-Hulbert with $\rho_\infty=1$, whose errors do not depend on
the start, and the error at $n=640$ is 603 (Chung-Hulbert, $\rho_\infty=0$) to
39513 (Newmark) times the error with the default start. With the default start
and the variable step sizes, the observed orders are 2.000 to 2.002. A step
retried after a solver failure gives the same result as the reduced step to the
last bit for all parameter sets.

The violation of the constraint is multiplied by $-\alpha_f/(1-\alpha_f)$ to
$10^{-12}$ in each step, i.e. by $0$, $-1$, $-0.8$, $-0.5$, $0$, $-0.125$ and
$0$ for Newmark, Chung-Hulbert with $\rho_\infty=1$, $0.8$, $0.5$, $0$, HHT and
WBZ. With the default start and $u(0)=(0.5,0)$, the constraint holds to
$10^{-12}$ in all 20 steps. The Newmark parameter set and `NewmarkBeta` agree to $1.2\cdot10^{-14}$.

**Reproducing the figures**

After building the `test_timestepmethods`,
`test_timestepmethods_stabilityregions`, `test_timestepmethods_symplectic`,
`test_multistepmethods`, `test_multistepmethods_stability` and
`test_generalizedalpha` executables, run:

```bash
python3 run_benchmarks.py <path-to-build-dir>
```

This runs all six executables to (re-)generate the JSON data and produces
the eight figures above via `plot_timestepmethods_convergence.py`,
`plot_timestepmethods_stabilityregions.py`,
`plot_timestepmethods_symplectic.py`, `plot_multistepmethods_convergence.py`,
`plot_multistepmethods_stability.py` and `plot_generalizedalpha.py`. If run
from within `<build-dir>/test/experimental/timestepping`, the build directory
argument can be omitted. Pass `--copy-to-doc` to additionally copy the figures
into `doc/doxygen/images/` (used to update the images shown on this page).
