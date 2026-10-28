# Benchmark Wrinkling of a Thin Sheet {#benchmark-sheet-wrinkling}

## Sheared Kapton Membrane (Wong and Pellegrino)

**Problem Description**

A thin sheet under in-plane load carries almost no compression: instead of shortening,
it buckles out of plane into a field of fine wrinkles whose wavelength is set by the
balance between the bending stiffness, which penalises short wavelengths, and the
tension along the wrinkles, which penalises long ones. This is the regime the
Föppl-von Kármán model (@ref FoepplVonKarmanPlate) is made for: the deflection is of
the order of the thickness, the rotations stay small, and the whole effect comes from
the quadratic slope term in the membrane strain.

The reference experiment is the sheared rectangular Kapton membrane of
Wong and Pellegrino \cite WongPellegrino2006. A rectangular foil is held along its two
long edges by rigid blocks; one block is displaced along the sheet by $\delta$ while the
distance between the blocks is kept fixed, and the two short edges are unrestrained.
Away from the short edges the state is one of simple shear, so the principal stresses
are $\pm G\gamma$ at $\pm 45^\circ$ with the shear strain $\gamma = \delta/W$. The
compressive principal stress wrinkles the sheet, and the wrinkles run at $45^\circ$
along the tensile direction.

**Parameters**

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Length | $L$ | 380 | mm |
| Height between the supports | $W$ | 128 | mm |
| Thickness | $t$ | 0.025 | mm |
| Young's modulus (Kapton HN) | $E$ | 3500 | N/mm² |
| Poisson ratio | $\nu$ | 0.31 | - |
| Shear displacement | $\delta$ | 0.5, 1, 2, 3 | mm |

**Reference Data**

Wong and Pellegrino scan the out-of-plane displacement along the mid-line
$y = 64$ mm and count each crest as one wrinkle:

| $\delta$ (mm) | wrinkles counted along the mid-line |
|---------------|-------------------------------------|
| 0.5 | 15 |
| 1.0 | 17 |
| 2.0 | 18 |
| 3.0 | 19 |

These are the counts of the loading branch: on unloading the sheet keeps the crests it
has, and at $\delta = 0.5$ mm it shows 19 of them where it showed 15 while loading. The
wrinkles are inclined at $45^\circ$ in the central region, with two "fan" regions near
the free short edges whose extent does not change with $\delta$. Their Figure 8 reports,
for the 0.025 mm foil, the average crest spacing along $x$ in the central region, which
the paper calls the wavelength $2\lambda$, and the average amplitude, read off the figure
here to about 0.3 mm and 0.005 mm:

| $\delta$ (mm) | 0.25 | 0.35 | 0.45 | 0.75 | 1.05 | 1.5 | 2.25 | 2.85 |
|---|---|---|---|---|---|---|---|---|
| crest spacing along $x$ (mm) | 32.3 | 29.2 | 27.0 | 25.0 | 24.0 | 22.0 | 21.1 | 20.5 |

| $\delta$ (mm) | 0.25 | 0.42 | 0.75 | 1.35 | 1.65 | 2.25 | 2.85 |
|---|---|---|---|---|---|---|---|
| amplitude (mm) | 0.150 | 0.200 | 0.250 | 0.297 | 0.330 | 0.365 | 0.382 |

The spacing along $x$ is the wavelength across the crests divided by $\cos 45^\circ$. Its
exponent in $\delta$, fitted over the measured range, is $-0.18$, and $-0.15$ above 0.45 mm: the
sheet keeps its crest count over a range of displacements and adds one in a jump.

**Two closed forms the wrinkled state must satisfy**

The bending-tension balance $\lambda = 2\pi (B/K)^{1/4}$ of Cerda and Mahadevan
\cite CerdaMahadevan2003, with $B = Et^3/(12(1-\nu^2))$ and the effective foundation
stiffness $K = N\pi^2/\ell^2$ of the tension $N$ acting along wrinkles of length
$\ell$, gives for the sheared sheet, where the wrinkles run diagonally so that
$\ell = \sqrt{2}W$ and $N = G\gamma t$,

$$\lambda = 2\pi\left(\frac{B\ell^2}{N\pi^2}\right)^{1/4}.$$

At $\delta = 0.5$ mm this is 21 mm. The length $\ell$ in it is an estimate, so the
balance fixes the exponent $\lambda \propto \delta^{-1/4}$ and the order of magnitude of
the prefactor, not its value.

The shear resultant is sharper. Once the wrinkles have relieved the compressive principal
stress, the sheet carries the shear as a uniaxial tension along the crests,
$N_1 = Et\,\gamma/2$, whose shear component is $N_{xy} = N_1/2$: the flat-sheet value
$G\gamma t$ times $(1+\nu)/2$, that is 0.655 for Kapton. A flat sheet, or one whose
wrinkles have not developed, carries the full $G\gamma t$.

With a strain $\varepsilon_0$ locked in across the sheet before the shear is applied, the
flat state is $\varepsilon_{xx} = 0$, $\varepsilon_{yy} = \varepsilon_0$,
$\varepsilon_{xy} = \gamma/2$, whose principal strains are $(\varepsilon_0 \pm R)/2$ with
$R = \sqrt{\varepsilon_0^2 + \gamma^2}$. The crests run along the tensile principal
direction, so all three closed forms generalise:

$$\tan 2\theta = -\frac{\gamma}{\varepsilon_0},\quad 2\theta \in (0,\pi),\qquad
\frac{N_{xy}}{G\gamma t} = \frac{1+\nu}{2}\left(1 + \frac{\varepsilon_0}{R}\right),\qquad
\varepsilon_c = \frac{(1-\nu)R - (1+\nu)\varepsilon_0}{2},$$

with $\theta$ the angle of the crests to $x$, $45^\circ$ at $\varepsilon_0 = 0$, and
$N = N_1 = Et(\varepsilon_0 + R)/2$ the tension along them in the balance above. At
$\varepsilon_0/\gamma = 0.1$ the crests turn by $2.9^\circ$, the shear resultant rises by
10 % and the strain the wrinkles absorb falls by 19 %, three signatures of the same number.
Tensioned this way and not sheared, the sheet carries the uniform strip solution
$N_{yy} = C\varepsilon_0$, $N_{xx} = \nu C\varepsilon_0$, $N_{xy} = 0$ away from the free
short edges, where $N_{xx}$ has to vanish. `Problem.PreStrain` sets $\varepsilon_0$ in the
scaled units, $(W/t)^2$ times the physical strain, as a prescribed $y$-displacement of the
moving edge.

**Scaled variables**

In units of $(W, t, t^2/W)$ for $(\mathbf{x}, w, \mathbf{u})$ the Föppl-von Kármán
system has no parameter but the Poisson ratio, the aspect ratio and the load, and every
field is of order one, which the numerically differentiated Jacobian needs. The input is
therefore $W = t = 1$, $E = 12(1-\nu^2)$ (so $D = 1$ and $C = 12$), $L/W = 2.96875$, and
the shear displacement scaled by $W/t^2 = 2.048\times10^5$ per mm, so that
$\delta = 0.5$ mm is $1.024\times10^5$. The scaled load is the distance from the buckling
threshold of the flat sheet, which is of order ten: the experiment sits four orders of
magnitude above it.

**Setup**

The domain is a structured triangulation of the rectangle. The Box scheme
(@ref BoxDiscretization) is used for the deformation sub-problem $(\varphi, w, \psi)$
and for the in-plane sub-problem $\mathbf{u}$, and the PQ1-Bubble scheme
(@ref PQ1BubbleDiscretization) for the rotations $\boldsymbol{\theta}$. The grid is a
simplex grid: with a single bubble the $\boldsymbol{\theta}/\psi$ pair is the MINI
element, which is stable on triangles but not on cubes.

The edges held by the rig are clamped, in plane and out of plane. The two short edges
are free: traction free in plane, and out of plane the natural condition of the mixed
formulation on all three deformation equations together with the traction
$\mathbf{T}\mathbf{n} = -\varphi\mathbf{n}$ on the rotations (see
@ref KirchhoffLovePlate). With a free edge present the shear curl potential is fixed by
that traction and carries no interior constraint. The two clamps are separate pieces of the
boundary joined by the free edges, so $\varphi = 0$ on the first clamp is the gauge while
$\varphi$ on the second is one unknown constant, imposed by the hole gauge of the model
(`Problem.ClampGauge = floating`, the default). Prescribing $\varphi = 0$ on both
(`fixed`) drops the net flux of $\nabla w - \boldsymbol{\theta}$ through the second clamp,
and both clamps may then release the sheet at the same uniform slope.

**Solution procedure**

The flat sheet is an exact solution at every load and, at these loads, a stationary point
of very high index: the linearised bending operator on the sheared state has about half
its eigenvalues negative. A Newton solve started anywhere near it seeks that saddle and
diverges, and loading up from the flat state chases the most unstable wavelength of the
flat sheet, which is an order of magnitude finer than the one the developed pattern
selects. The run therefore starts from the tension-field state
(`Problem.SeedTensionField`): a $45^\circ$ wrinkle field at the wavelength of the balance,
or at the wavelength `Problem.SeedWavelength`, whose amplitude
$A = 2\sqrt{(1-\nu)\gamma/2}/k$ absorbs the part of the compressive principal strain
$\gamma/2$ that the Poisson contraction of the tension along the crests does not supply,
together with the in-plane displacement $u_n = -(A^2k/8)\sin(2\mathbf{k}\cdot\mathbf{x})$ that cancels the
part of the wrinkle strain oscillating at twice the wavenumber. Without the second the
seed carries more spurious membrane stress than the applied load. The rotations start
from the slope of the seeded deflection and the shear gradient potential from
$\varphi = Dk^2 w$, which is what the rotation balance $\nabla\varphi = -D\nabla\Delta w$
gives for a wave of wavenumber $k$: the transverse balance then starts with the bending
force of the seed to absorb instead of its membrane traction, larger by $N/(Dk^2)$. On
216 x 72 cells at $\delta = 0.5$ mm the first pseudo-time step converges only below
$10^{-10}$ without the potential and every attempt to grow it failed; with it the step
reached $3\times10^{-6}$ after thirty steps. The seed's envelope across the sheet is a
boundary layer $\tanh^2(y/\ell)\tanh^2((W-y)/\ell)$ with $\ell = \lambda/\pi$, which has
zero slope at the clamped edges.

From the seed the state is relaxed in pseudo-time: $c\,\partial_\tau w = -\partial E/\partial w$
integrated implicitly descends the energy and can only come to rest at a stable
equilibrium. The relaxation stops once the deflection moves by less than
`Relaxation.SteadyStateTolerance` thicknesses in a step of the largest size allowed. The
first step from the seed has to be shorter than the relaxation time $1/(Dk^4)$ of the
stiffest wrinkle modes, about $10^{-7}$ here, or it is a static solve from an
inconsistent state: Newton diverges, and the direct solver still has to factorise the
tangent of the diverged iterate. `Relaxation.InitialTimeStepSize` is
$10^{-9}$, the step controller grows the step by about 1.6 per converged step and reaches
the cap `Relaxation.MaxTimeStepSize` of 0.05 after some forty steps, and a later load of a
continuation starts from the relaxed state at `Relaxation.ContinuationTimeStepSize`,
$10^{-5}$. A static solve is not attempted from the relaxed state: the pattern can slide along its own
normal at almost no cost, so the undamped tangent is nearly singular there and Newton
wanders off (the first iteration raises the residual by a factor of 400).

Newton runs without a line search. The residual norm is dominated by the in-plane
equations, whose linearisation error after a step of the deflection is quadratic in that
step and larger than the small residual each pseudo-time step starts from, so a line
search that insists on a monotone residual rejects full steps for ten iterations before
the quadratic convergence sets in; on the regression case a step takes 13 iterations with
the line search and 8 without, and the two Newton paths relax onto different metastable
states, the one without the line search 0.06 % lower in energy.

**Results**

The run writes `<name>_wrinkles.dat` with the amplitude over the central half of the
mid-line, the crest count, the crest spacing, the wrinkle orientation from the structure
tensor of the slope, the wavelength, the central shear resultant against the flat and the
tension-field values, the number of relaxation steps and the bending and membrane strain
energy, which is what tells two metastable states apart; `<name>_profile_*.dat` holds the
scan line and the VTK output the deflection and the stress resultants. During the
relaxation every step prints the pseudo-time step, the largest deflection, the change of
the deflection, the central amplitude and the strain energy.

The ctest runs $\delta = 20\,\mu$m ($4096$ scaled) on 48 x 16 cells, reached in 88
relaxation steps:

| quantity | measured | reference |
|---|---|---|
| central shear resultant | 1.020 $\times$ tension field | $(1+\nu)/2 = 0.655$ of the flat value |
| wrinkle normal | $-45.5^\circ$ | $-45^\circ$ |
| wavelength | 0.347 $W$ | 0.369 $W$ from the balance |
| amplitude over the central mid-line | 5.5 $t$ | developed, $A/t \gg 1$ |

and asserts the shear resultant to 5 %, the orientation to $10^\circ$, the wavelength
to 30 % and an amplitude above two thicknesses. The shear resultant is the sharp one of
the four: it distinguishes the relaxed tension-field state from the flat sheet by a factor
$(1+\nu)/2$, and the measured value is 2 % above it. On 60 x 20 cells the same load gives
the normal at $-45.3^\circ$, $\lambda = 0.325\,W$ and the shear resultant at 1.018 of the
tension-field value; there the free-edge corners, which curl up to 16 thicknesses, keep
moving in discrete coarsening events after the central pattern has settled, and the
per-step stopping rule is not met within 150 steps.

Two further configurations gate the strain locked in across the sheet, one signature of
$\varepsilon_0$ each; the third, the strain the wrinkles absorb, sets the amplitude of the
seed and not of the relaxed state, whose mid-line scan crosses the crests where the envelope
along them peaks and so runs a fifth above the uniform value, and is therefore not gated. `test_fvk_wrinkling_prestrain_flat` tensions the sheet by
$\varepsilon_0 = 4\times10^{-4}$ and does not shear it, 96 x 32 cells, four Newton
iterations: over the central half the resultant across the sheet is 0.985 of
$C\varepsilon_0$, asserted to 3 %, the one along it 0.84 of $\nu C\varepsilon_0$ and the
shear resultant $6\times10^{-6}$ of $C\varepsilon_0$. The shortfall in $N_{xx}$ is not an
error: the free short edges relieve it over a length of the order of $W$, which reaches the
middle of a sheet only $2.97\,W$ long, and that is why the gate is on $N_{yy}$.
`test_fvk_wrinkling_prestrain` repeats the regression with $\varepsilon_0/\gamma = 0.1$:
the wrinkle normal is $-42.8^\circ$ against the $-42.1^\circ$ of the closed form and
the shear resultant 1.017 of the tension-field value, asserted to $2^\circ$ and 5 %. The
same state without the locked-in strain sits at $-45.5^\circ$, so the measured normal moves
by $2.7^\circ$ between the two, more than the tolerance: a strain that did not reach the
wrinkled state would fail the gate.

**The experimental displacement of 0.5 mm: one displacement, several states**

At $\delta = 0.5$ mm the sheet has a metastable state for every crest count in a band, and
a relaxation ends in the state nearest to where it starts: from the balance seed on
108 x 36 cells the crest spacing after 100 steps was that of the seed to 2 %. The states
are therefore compared by their strain energy. `wrinkling_benchmark.py --seed-crests`
seeds $n$ crests along the mid-line; on 144 x 48 cells (2.7 mm) with 300 relaxation steps
this gives, with the number of maxima above a quarter of the central amplitude, the mean
crest spacing along $x$ over the central half, the central amplitude, the wavelength across
the crests and the strain energy relative to the lowest:

| seeded | counted | spacing (mm) | $A$ (mm) | $\lambda$ (mm) | $E/E_{\min} - 1$ |
|---|---|---|---|---|---|
| 11 | 12 | 28.1 | 0.288 | 20.2 | $1.5\times10^{-4}$ |
| 13 | 12 | 28.1 | 0.282 | 20.2 | $1.2\times10^{-4}$ |
| 15 | 12 | 25.6 | 0.260 | 18.4 | 0 |
| 17 | 14 | 22.6 | 0.231 | 16.2 | $1.9\times10^{-4}$ |
| 19 | 15 | 20.5 | 0.236 | 14.7 | $5.5\times10^{-4}$ |

The 15 seed relaxes to twelve crests between $x = 50$ and 330 mm at 25.6 mm with the
corners curled up to 0.8 mm. It is the lowest in energy of the five, its spacing is 4 %
below the 27 mm of Figure 8 at $\delta = 0.45$ mm, and 380 mm at that spacing hold the 15
crests Figure 4 counts. The 11 and 13 seeds relax to the same 28.1 mm state, which is what
the balance seed produces, $1.2\times10^{-4}$ above the lowest with the corners in two
slightly different configurations; 14 and 15 central crests are metastable at
$1.9\times10^{-4}$ and $5.5\times10^{-4}$ above, with the corners curled the other way,
the last at the 20.5 mm spacing the experiment shows at 2.85 mm and after unloading.
Which basin a seed ends in depends on the seed's details; the states and their energies
do not. The band of
metastable counts is what makes the measured count path dependent. On 216 x 72
cells the 13, 15 and 17 seeds relax to 28.7, 24.9 and 22.9 mm with the energies
$2.2\times10^{-4}$, 0 and $1.1\times10^{-4}$ above the lowest, so the ordering and the
spacing of the lowest state, 7 % below the measured, hold on the finer mesh, which lowers
every energy by $2.6\times10^{-3}$. The energies
settle to $10^{-6}$ per step; only the 13 seed met the per-step rule on the deflection
within 300 steps, the others still drifted by 0.02-0.05 thicknesses per step at the
corners. The mid-line amplitude 0.260 mm sits inside the 0.25-0.3 mm of the central
peaks in the scan of Figure 4(a) and 22 % above the average of Figure 8; the uniform
tension-field amplitude $(\lambda/\pi)\sqrt{(1-\nu)\gamma/2}$ at the computed wavelength
is 0.214 mm, and the scan line at $y = W/2$ crosses the crests where a half-sine envelope
along them peaks. All five states have the wrinkle normal at $-44.1^\circ$ and the central
shear resultant at 0.974 of the tension-field value. The fans are multistable as well:
in most of these states the corners curl monotonically over the outer 50 mm to 0.8-0.9 mm,
whereas the 13-crest state on 216 x 72 cells has flat ends and two peaks in each fan, at 12
and 49 mm from one edge and 331 and 368 mm from the other, which is the structure Figure 4
shows, with the measured ends at $-0.6$ to $-0.8$ mm.

The rig holds the sheet horizontally, so its weight, 0.35 Pa for the 25 $\mu$m foil, acts
out of plane; in the scaled units that is a transverse load `Problem.Gravity` of
$qW^4/(Dt) = 7.4\times10^5$. Applied with either sign relative to the 15-crest seed it
leaves the central pattern unchanged, spacing 25.3-25.5 mm, normal $-44.1^\circ$, energy
within $10^{-5}$, and acts on the slack corners alone: both free edges go 1.5-1.6 mm along
the load instead of curling by 0.8, and one hump of 0.4 mm appears 34 mm from each edge.
That is the sign pattern of the measured ends, at twice their excursion and with one peak
instead of two, so the weight alone does not produce the measured fans.

**The strain locked in before the shear**

Before the fixed edge is clamped the rig tensions the sheet in $y$ with counterweights, a
"small initial stress" the paper does not quantify, so the sheet carries a strain
$\varepsilon_0$ across its width under every shear displacement. `Problem.PreStrain` sets
it, in the scaled units $(W/t)^2\varepsilon_0$, as a prescribed $y$-displacement of the moving
edge; the seed absorbs $(1-\nu)\gamma/2 - (1+\nu)\varepsilon_0/2$ and the tension-field
resultant carries the factor $1 + \varepsilon_0/\gamma$. Such a strain rotates the tensile
principal direction and with it the crests away from $45^\circ$ towards $y$ by an angle that
falls with $\gamma$, which the photographs of the experiment can check: the crest angle in
their Figure 3, read by a gradient-orientation histogram over the central half of each
panel and corrected for the stretch of the panel on the page, is $48.3^\circ$, $46.4^\circ$,
$46.3^\circ$ and $46.2^\circ$ at 0.5, 1, 2 and 3 mm, the first to about $\pm2^\circ$ since
its wrinkles are faint. On 144 x 48 cells from the 15-crest seed at 0.5 mm:

| $\varepsilon_0$ | crest angle | spacing (mm) | $A$ (mm) | ends of the scan (mm) | fan maxima at (mm) |
|---|---|---|---|---|---|
| 0 | $45.9^\circ$ | 25.6 | 0.260 | +0.81 / +0.64 | none |
| $2\times10^{-4}$ | $47.1^\circ$ | 24.5 | 0.256 | +0.72 / $-0.75$ | 48, 348 |
| $4\times10^{-4}$ | $48.3^\circ$ | 24.9 | 0.244 | +0.23 / $-0.22$ | 8, 45, 348 |
| $8\times10^{-4}$ | $50.7^\circ$ | 23.4 | 0.210 | $-0.04$ / +0.05 | 26, 340, 369 |
| photographs and scans | $48.3^\circ$ | 26.7 | 0.21 | $-0.6$ to $-0.8$ | two per fan |

$\varepsilon_0 = 4\times10^{-4}$, a resultant of 0.039 N/mm across the sheet against the
0.13 N/mm of shear at 0.5 mm, gives the photographed angle at 0.5 mm, and from their own
balance seeds $47.1^\circ$ at 1 mm and $46.6^\circ$ at 2 mm against $46.4^\circ$ and
$46.3^\circ$ photographed, where the model without it stays at $45.9^\circ$; the spacing
is then 23.7 and 20.5 mm against the measured 24.0 and 21.3, and the amplitude 0.325 and
0.425 mm against 0.29 and 0.36. In the fans the tension across the sheet stops the corners
curling: the outermost crests become the largest, 0.4 to 0.5 mm at 0.5 mm and 0.6 to 0.9
mm at 1 and 2 mm, and at 1 mm each fan carries two peaks, at 16 and 48 mm from one edge
and 332 and 364 mm from the other, which is the structure of the measured scans. What
this strain does not give is the measured excursion of $-0.6$ to $-0.8$ mm at both ends of
the scan; the weight does. With both, $\varepsilon_0 = 4\times10^{-4}$ and the weight of
the horizontal sheet, the scan at 0.5 mm has both ends at $-0.72$ and $-0.69$ mm, two peaks
of 0.39 and 0.42 mm in each fan at 29 and 45 mm from the edge, wrinkles of 0.2 to 0.25 mm
between them at a spacing of 24.1 mm and the crests at $48.3^\circ$, with the central
amplitude 0.229 mm: every feature of the measured scan, the spacing 10 % finer. The seed
of the opposite sign relative to the weight gives the mirror image with both ends up. At
1 mm the same pair gives the ends at $-0.40$ and $-0.68$ mm, one peak of 0.6 mm per fan,
a spacing of 23.4 mm against the measured 24.0, the crests at $47.1^\circ$ and the
central amplitude 0.329 mm against 0.29.

**Loading and unloading**

The experiment loads the sheet from zero and counts more crests at every displacement, then
keeps them on unloading. `wrinkling_benchmark.py --path` runs the displacements in sequence,
each relaxed from the state of the previous one with the in-plane field beyond the locked-in
strain rescaled to the new displacement, here 200
steps per displacement on 144 x 48 cells from the balance seed at the first; the measured
values are Figure 8 interpolated in $\log\delta$, and 0.38 mm at 3 mm is its last point at
2.85 mm:

| $\delta$ (mm) | path | counted | spacing (mm) | measured | $A$ (mm) | measured | $\lambda$ (mm) | balance |
|---|---|---|---|---|---|---|---|---|
| 0.1 | loading | 8 | 36.9 | | 0.178 | | 26.3 | 31.6 |
| 0.25 | loading | 10 | 33.2 | 32.3 | 0.238 | 0.150 | 23.8 | 25.1 |
| 0.5 | loading | 12 | 28.1 | 26.7 | 0.284 | 0.21 | 20.2 | 21.1 |
| 1 | loading | 13 | 25.6 | 24.0 | 0.376 | 0.29 | 18.4 | 17.8 |
| 2 | loading | 16 | 19.5 | 21.3 | 0.415 | 0.36 | 14.0 | 15.0 |
| 3 | loading | 17 | 19.5 | 20.5 | 0.502 | 0.38 | 14.0 | 13.5 |
| 2 | unloading | 16 | 18.2 | 21.3 | 0.384 | 0.36 | 13.1 | 15.0 |
| 1 | unloading | 16 | 18.2 | 24.0 | 0.280 | 0.29 | 13.1 | 17.8 |
| 0.5 | unloading | 14 | 22.6 | | 0.233 | 0.21 | 16.2 | 21.1 |

On loading the spacing is 5 and 7 % above the measured at 0.5 and 1 mm and 8 and 5 % below
at 2 and 3 mm, the count jumping from 13 to 16 between 1 and 2 mm where the experiment goes
from 17 to 18; over 0.25 to 3 mm the spacing falls with the exponent $-0.23$ against the
measured $-0.18$ and the $-0.25$ of the balance. Unloading keeps the pattern of 3 mm down to
1 mm, 16 crests at 18.2 mm, as the experiment keeps its 19 crests down to 0.5 mm, and at
0.5 mm settles in the 14-crest state of the table above, so between loading and unloading
the spacing at 0.5 mm spans 28.1 to 22.6 mm against the 15 and 19 crests of the experiment.
At 1 and 2 mm the unloading states lie $8\times10^{-5}$ and $5\times10^{-5}$ below the
loading states in energy: the finer pattern is the lower one and loading follows a
metastable branch, as the measured counts do. The amplitude on the mid-line runs 5 to 30 %
above the averages of Figure 8 on loading and within 10 % on unloading; the root-mean-square
amplitude of the scan line is 5 to 10 % below the half peak-to-peak, and the tension-field
amplitude at the computed wavelength, $(\lambda/\pi)\sqrt{(1-\nu)\gamma/2}$, is 0.21, 0.30,
0.32 and 0.40 mm at 0.5, 1, 2 and 3 mm, within 10 % of Figure 8, so the amplitude
comparison mostly repeats the one of the spacing. Relaxed from their own balance seeds on
216 x 72 cells, the displacements 1, 2 and 3 mm give 23.1, 19.7 and 18.6 mm, 4 to 9 %
below the measured, with 14, 16 and 17 crests counted.

The `Stretch` loading option runs the other classical setup instead, the sheet pulled
along its length between two grips, in which the clamped grips prevent the lateral
Poisson contraction and put the centre of the sheet into transverse compression
\cite CerdaMahadevan2003. A linear-material Föppl-von Kármán model gives a monotone
pitchfork for that problem and therefore cannot reproduce the flattening of the wrinkles
at large strain that hyperelastic models and the experiments show \cite LiHealey2016.
