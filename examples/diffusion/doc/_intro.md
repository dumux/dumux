# Diffusion equation

This example is implemented in three files, `model.hh`, `params.input`, `main.cc`.
The `model.hh` header implements a simple self-contained diffusion equation model
to be used with a control-volume finite element discretization like the Box method.

__The main points illustrated in this example are__

* Setting up a simple model
* Random initial solution in parallel
* Solving a time-dependent diffusion problem
* Use the linear PDE solver
* Reuse the Jacobian for all time steps

[TOC]

## Equation and problem description

The scalar diffusion equation on a domain $\Omega \subset \mathbb{R}^2$
with boundary $\partial\Omega = \Gamma_D \cup \Gamma_N$ composed of Dirichlet and Neumann boundaries
reads

```math
\begin{aligned}
\frac{\partial c}{\partial t} - \nabla \cdot (D \nabla c) &= 0 \quad \mathrm{in}\; \Omega \\
                                                        c &= c_D  \quad \mathrm{on}\; \Gamma_D \\
                         -D \nabla c \cdot \boldsymbol{n} &= g_N    \quad \mathrm{on}\; \Gamma_N \\
\end{aligned}
```

with the unknown concentration $c(x,t)$, diffusion coefficient $D$, Dirichlet data $c_D$, Neumann data $g_N$,
and $\boldsymbol{n}$ denoting the outward-oriented unit normal vector on $\partial\Omega$.
In this example we will use homogenenous Neumann boundary conditions
on all of $\partial\Omega$ such that $\Gamma_D = \emptyset$ and $g_N = 0$.

We partition $\Omega$ into control volumes $B$ (sometimes called "boxes"). Integrating our equation over a control volume $B$
and applying the divergence theorem yields

```math
\begin{aligned}
\int_B \frac{\partial c}{\partial t} \mathrm{d}V - \int_{\partial B} D (\nabla c \cdot \boldsymbol{n}) \mathrm{d}A &= 0.
\end{aligned}
```

We discretize in time using an implicit Euler method

```math
\begin{aligned}
\int_B \frac{c^{n+1}-c^{n}}{\Delta t} \mathrm{d}V - \int_{\partial B} D (\nabla c^{n+1} \cdot \boldsymbol{n}) \mathrm{d}A &= 0.
\end{aligned}
```
The initial data is given by $c^0 := c(x, 0)$. In our example, $c^0$ is a random field, where the control volume
averages $c^0_{h,B}$ are sampled from a uniform distribution, $c^0_{h,B} \sim U(0,1)$.

For the Box method, the control volumes are centered at grid vertices. We partition the control volume boundary
into sub control volume faces $\sigma$ such that each $\sigma$ belongs to exactly one grid element. On each element $K$,
the discrete solution is represented by a linear basis functions as

```math
c_{h,B}(x) = \sum\limits_{K \in \mathcal{B}_K} c_{h,B} N_K(x)
```

where $\mathcal{B}_K$ is the set of control volumes intersecting element $K$ and $N_B$ are the basis functions associated
with the vertex corresponding to the center of control volume $B$. Insertion of this ansatz yields the discrete equation

```math
\begin{aligned}
\vert B \vert \frac{c_B^{n+1}-c_B^{n}}{\Delta t} - \sum_{\sigma \in \Sigma_{BK}} \left[ D \sum_{B \in \mathcal{B}_K} c^{n+1}_{h,B} \nabla N_B \cdot \boldsymbol{n}_{B,\sigma} \vert \sigma \vert \right] &= 0,
\end{aligned}
```

where $\Sigma_{BK}$ is the set of sub control volume faces which belong to $B$ and $K$, and $\boldsymbol{n}_{B,\sigma}$ is the outward-oriented
unit normal vector on $\partial B$.

### Model parameters

* $\Omega = [0,1]\times[0,1]$ (unit square)
* $T = [0,5]$ with time step size $\Delta t =  0.1$
* $100 \times 100$ cells (structured Cartesian grid)
* $D = 0.0001$
* $c^0_{h,B} \sim U(0,1)$

### Simulation result

The simulation result will look something like this.

<figure><img src="img/diffusion.gif" alt="diffusion"/></figure>

### Running the example in parallel

By default Dumux will try to speed up the assembly by using shared memory parallelism if a suitable
backend has been found on your system (one of TBB, OpenMP, Kokkos, C++ parallel algorithms).
You can limit the number of threads by prepending your executable with `DUMUX_NUM_THREADS=<number>`.
If you also want to use distributed memory parallelism with MPI (works better for solvers at the moment),
run the executable with your MPI environment. Each MPI process will use multi-threading if
`DUMUX_NUM_THREADS` is larger than $1$.

Running the example with four MPI processes (distribution memory parallelism)
each with two threads (shared memory parallelism):

```sh
DUMUX_NUM_THREADS=2 mpirun -np 4 ./example_diffusion
```

You can set the parameter `Grid.Overlap` to some non-zero integer in `params.input`
to turn the domain decomposition into an overlapping decomposition where
`Grid.Overlap` specifies the number of grid cells in the overlap between processes.
This can help to increase the convergence speed of the linear solver.

# Implementation

For the code implementation see Part I and Part II.
