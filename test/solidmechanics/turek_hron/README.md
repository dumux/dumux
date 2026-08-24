# CSM Benchmark (Turek & Hron) {#benchmark-turekhron-csm}

This benchmark test case implements the Computational Solid Mechanics (CSM) benchmark as proposed by Turek and Hron @cite TurekHron2006. It evaluates the performance of the hyperelastic model for large deformations and time-dependent loading.

## Problem description

The setup consists of a fixed rigid cylinder with an elastic beam (flag) attached to its rear. The beam is subject to a constant volume force (gravity), which induces oscillation. The benchmark tests the ability of the numerical scheme to handle large displacements and the accuracy of the time integration.

## Mathematical model

We use the hyperelastic model @ref Hyperelastic. This requires implementing the first Piola stress tensor $\mathbf{P}$ as a function of the deformation gradient $\mathbf{F} = \mathbf{I} + \nabla \mathbf{u}$. \
We use a **Saint-Venant Kirchhoff material law** similar to the one described in  @cite TurekHron2006, which defines the first Piola-Kirchhoff stress tensor as:
$$
\mathbf{P}(\mathbf{F}) = \mathbf{F} \mathbf{S}
$$
with
* **Lagragian-Green strain tensor**: $\quad\mathbf{E} = 0.5(\mathbf{F}^T\mathbf{F} - \mathbf{I}) = 0.5(\mathbf{C} - \mathbf{I})$
* **Second Piola-Kirchhoff stress tensor**: $\quad\mathbf{S}(\mathbf{E}) = \lambda \text{tr}(\mathbf{E})\mathbf{I} + 2\mu\mathbf{E}$,

where $\lambda$ is the first Lamé parameter and $\mu$ is the shear modulus.

### Numerical Discretization
| Implementation | Spatial Discretization | Type | Order | Temporal Discretization | Type | Order |
| :--- | :--- | :---: | :---: | :--- | :--- | :---: |
| Dumux | Box @ref BoxDiscretization | CVFE | 1 | Newmark-$\beta$ | Implicit | 2 |
| Reference (Featflow) | Galerkin FEM  (Q2/P1) | FEM | 2 | Crank-Nicolson | Implicit | 2 |

Due to the difference in interpolation order (linear vs. quadratic), DuMuX typically requires a higher mesh refinement level to match the peak amplitudes of the reference solution.

## Setup & Parameters

The Turek and Hron suite outlines three specific structural sub-benchmarks. While this implementation focuses primarily on **CSM3** (the instationary time-dependent case), the domain configuration and parameters across the suite are summarized below for complete context.

### Domain Geometry and Boundary Segments

![CSM Geometry and Boundary Conditions](solid_turekhron_domain.svg)

The domain isolates the flexible structural appendage from the benchmark layout using key geometric references:

| Parameter | Symbol | Value | Unit |
| :--- | :---: | :---: | :---: |
| Total length of the hyperelastic beam | $L$ | $0.35$ | $\text{m}$ |
| Total height thickness of the hyperelastic beam | $H$ | $0.02$ | $\text{m}$ |
| Initial coordinates of monitoring control point $A$ (beam tip center) | $A$ | $(0.6, 0.2)$ | $\text{m}$ |
| Coordinates of reference boundary point $B$ (upstream cylinder hull face) | $B$ | $(0.15, 0.2)$ | $\text{m}$ |
| Coordinates of reference point $C$ (rigid cylinder base center) | $C$ | $(0.2, 0.2)$ | $\text{m}$ |

### Boundary and Initial Conditions

The physical behavior is constrained by prescribing conditions along the closure of the solid boundary $\partial\Omega = \Gamma_D \cup \Gamma_N$:

* **Dirichlet Boundary ($\Gamma_D$)**: The left vertical face of the elastic structure is clamped flat against the rigid cylinder profile. The displacement field is strictly fixed:
  $$\mathbf{u} = \mathbf{0} \quad \text{on } \Gamma_D \quad (x = 0.2\text{ m})$$

* **Flux Boundary ($\Gamma_N$)**: All remaining outer surfaces (top, bottom, and right end) exposed to environmental conditions are completely traction-free. The exact mathematical formulation defines a zero surface traction force vector using the first Piola-Kirchhoff stress tensor $\mathbf{P}$ and the outward unit normal vector $\mathbf{n}$:
  $$\mathbf{P} \cdot \mathbf{n} = \mathbf{0} \quad \text{on } \Gamma_N$$

* **Initial condition**: The domain begins in a fully relaxed, non-deformed state at time $t = 0$. For all steps $t > 0$, a constant spatial gravity vector ($g = 2\text{ m/s}^2$) acts downward along the negative $y$-axis across the solid density.

### Parameter settings for the CSM tests

| Parameter | Symbol | Value (CSM1) | Value (CSM2) | Value (CSM3) | Unit |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Solid Density | $\rho^s$ | $1.0 \times 10^3$ | $1.0 \times 10^3$ | $1.0 \times 10^3$ | $\text{kg/m}^3$ |
| Poisson's Ratio | $\nu^s$ | $0.4$ | $0.4$ | $0.4$ | $-$ |
| Shear Modulus | $\mu^s$ | $0.5 \times 10^6$ | $2.0 \times 10^6$ | $0.5 \times 10^6$ | $\text{kg/(m}\cdot\text{s}^2\text{)}$ |
| Young's Modulus | $E$ | $1.4 \times 10^6$ | $5.6 \times 10^6$ | $1.4 \times 10^6$ | $\text{Pa}$ |
| Gravity | $g$ | $2.0$ | $2.0$ | $2.0$ | $\text{m/s}^2$ |
| Simulation Type | $-$ | Steady State | Steady State | Instationary | $-$ |
And the corresponding Lamé parameters are computed as:
$$
\lambda = \frac{\nu^s E}{(1+\nu^s)(1-2\nu^s)} \quad \mu^s = \frac{E}{2(1+\nu^s)}
$$

## Data creation and comparison

The primary objective is to monitor the displacement fields $u_x$ and $u_y$ of control point A (over time).
To verify the accuracy of this numerical implementation, the values can be validated against the highly refined grid-independent reference results provided by the TU Dortmund Featflow benchmark suite [1].

### Steady State Tests (CSM1 & CSM2)
If utilizing this setup to evaluate steady-state accuracy before running time-dependent loops, compare your final deflection values against the highest mesh-refinement levels:

![CSM1 and CSM2 Validation](csm12_comparison.png)
The plots above compare DuMuX (Box scheme) against the Featflow reference data for stationary cases.

### Instationary Test (CSM3)
For the dynamic **CSM3** case, the beam exhibits large, persistent oscillations due to the sudden application of gravity on the elastic structure starting from an undeformed configuration.

![CSM3 Validation](csm3_comparison.png)

| $\Delta t$ | nel | ndofs | Source | Discretisation | ux of A [mm] (Mean $\pm$ Amp [Freq]) | uy of A [mm] (Mean $\pm$ Amp [Freq]) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.005 | 320 | 6468 | FeatFlow | FEM | -14.279 +/- 14.280 [1.0995] | -63.541 +/- 65.094 [1.0995] |
| 0.005 | 1280 | 25092 | FeatFlow | FEM | -14.299 +/- 14.299 [1.0995] | -63.594 +/- 65.154 [1.0995] |
| 0.005 | 5120 | 98820 | FeatFlow | FEM | -14.305 +/- 14.305 [1.0995] | -63.607 +/- 65.160 [1.0995] |
| 0.010 | 320 | 6468 | FeatFlow | FEM | -14.632 +/- 14.636 [1.0978] | -64.744 +/- 64.907 [1.0978] |
| 0.010 | 1280 | 25092 | FeatFlow | FEM | -14.645 +/- 14.650 [1.0978] | -64.765 +/- 64.946 [1.0978] |
| 0.010 | 1923 | 1087 | DuMuX | Box | -13.719 +/- 13.720 [1.2295] | -62.448 +/- 63.292 [1.1111] |
| 0.010 | 5120 | 98820 | FeatFlow | FEM | -14.645 +/- 14.650 [1.0978] | -64.766 +/- 64.948 [1.0978] |
| 0.010 | 10631 | 5586 | DuMuX | Box | -14.396 +/- 14.402 [1.2097] | -63.970 +/- 64.709 [1.0949] |
| 0.010 | 49710 | 25429 | DuMuX | Box | -14.590 +/- 14.594 [1.2097] | -64.381 +/- 65.109 [1.0949] |
| 0.010 | 218076 | 110216 | DuMuX | Box | -14.629 +/- 14.633 [1.2097] | -64.497 +/- 65.178 [1.0989] |
| 0.020 | 320 | 6468 | FeatFlow | FEM | -14.384 +/- 14.389 [1.0956] | -64.271 +/- 64.595 [1.0956] |
| 0.020 | 1280 | 25092 | FeatFlow | FEM | -14.402 +/- 14.406 [1.0956] | -64.352 +/- 64.679 [1.0956] |
| 0.020 | 1923 | 1087 | DuMuX | Box | -13.659 +/- 13.665 [1.2295] | -61.923 +/- 63.544 [1.1194] |
| 0.020 | 5120 | 98820 | FeatFlow | FEM | -14.404 +/- 14.408 [1.0956] | -64.371 +/- 64.695 [1.0956] |
| 0.020 | 10631 | 5586 | DuMuX | Box | -14.336 +/- 14.338 [1.2097] | -63.910 +/- 64.510 [1.1029] |
| 0.020 | 49710 | 25429 | DuMuX | Box | -14.516 +/- 14.521 [1.2097] | -64.315 +/- 64.933 [1.0870] |
| 0.020 | 218076 | 110216 | DuMuX | Box | -14.580 +/- 14.584 [1.1905] | -64.381 +/- 65.104 [1.0870] |

### Reproducing our simulation data

To reproduce the benchmark results, navigate to the specific test directory within your build folder and run the benchmark script.

```bash
# Example for CSM12
cd build-cmake/test/solidmechanics/turek_hron/csm12
make test_saint_venant_kirchhoff
python3 runBenchmarks.py
# Example for CSM3
cd build-cmake/test/solidmechanics/turek_hron/csm3
make test_saint_venant_kirchhoff_dynamic
python3 runBenchmarks.py
```

## References
[1] Featflow Benchmark Suite, TU Dortmund. *CSM Tests Reference Sheet*. Available online at: https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/fsi_benchmark/fsi_tests/fsi_csm_tests.html
