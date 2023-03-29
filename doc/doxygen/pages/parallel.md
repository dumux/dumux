# Running in parallel

In this section, we discuss how to run DuMux applications in parallel.
There is two types of parallelism implemented in DuMux which can be combined: shared memory parallelism (multi-threading)
and distributed memory parallelism.

[TOC]

## Distributed memory parallelism with MPI

Distributed memory parallelism with MPI (Message Passing Interface) is supported through DUNE. The main idea behind the MPI parallelization is the concept of domain decomposition. For parallel
simulations, the computational domain is split into subdomains and one process (rank) is used to solve the local problem of each subdomain.
During the global solution process, some data exchange between the ranks/subdomains is needed.
MPI is used to send data to other ranks and to receive data from other ranks.

Domain decomposition in DUNE is handled by the grid implementations.
The grid is partitioned and distributed onto several processes. The `grid.leafGridView()` and `grid.levelGridView()` methods return views on the processor-local partition of the grid. For example, `grid.leafGridView().size(0)` returns the number of elements available on the current processor.
Most grid implementations contain their own domain decomposition methods to split the computational domain into subdomains. Some grid implementations also support external tools like `METIS`, `ParMETIS`, `PTScotch` or `ZOLTAN` for partitioning. These tools may have to be installed separately. Consult the documentation of the grid implementation.

Linear algebra types such as matrices and vectors do not know that they are in a parallel environment. This means, for example, that the solution vector on each processor only contains the number of degrees of freedom available on that processor and its size is therefore smaller than the global number of degrees of freedom.
Communication of linear algebra type data is handled by parallel solvers.
Currently, most solvers in DuMux are also usable as parallel solvers.

The solvers based on dune-istl can be found in dumux/linear/istlsolvers.hh
and dumux/linear/istlsolverfactorybackend.hh. The latter makes the solver
configurable via parameters or the command line interface. To use the dune-istl based iterative solvers in parallel, you have to use the constructor receiving a grid view and a dof mapper, and the
linear solver traits derived from the grid geometry type, e.g.

```cpp
#include <memory>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

int main {

using GridGeometry = ...;
using Assembler = ...;
auto gridGeometry = std::make_shared<GridGeometry>(...);

...

// type of the Linear solver
using LinearSolver = AMGBiCGSTABIstlSolver<
    LinearSolverTraits<GridGeometry>>,
    LinearAlgebraTraitsFromAssembler<Assembler>
>;

// construct a parallel Linear solver
auto linearSolver = std::make_shared<LinearSolver>(
    gridGeometry->gridView(), gridGeometry->dofMapper()
);

...
}
```

In order for DuMux simulation to run in parallel, an MPI library (e.g. OpenMPI, MPICH or IntelMPI) implementation must be installed on the system.

Note that the parallel AMG preconditioner of dune-istl defaults
to an iterative SSOR coarse grid solver if no direct solver is found on your system. Unfortunately, the iterative solver has a very high and hard-coded tolerance as a termination criterion, which will not solve the coarse grid system with sufficient accuracy for typical problems in DuMux.
We therefore recommend to install one of the direct solver libraries supported by dune-istl. This is either UMFPack contained in SuiteSparse, or SuperLU, see also the section on [External Libraries](#external-libraries).

## Run a parallel MPI application

The starting procedure for parallel simulations depends on the chosen MPI library. Most MPI implementations use the `mpirun` command

```sh
   mpirun -np [n_cores] [executable_name]
```

where `-np [n_cores]` sets the number of cores that should be used for the computation.
On a cluster you usually have to use a queuing system (e.g. slurm) to submit a job.
Check with your cluster administrator how to run parallel applications on the cluster.

## Handling parallel results (VTK output)

For serial computations, DuMux produces single vtu-files as default output format.
During a simulation, one `*.vtu`/`*.vtp` file is written for every output step.
In the parallel case, one `*.vtu`/`*.vtp` file for each step and processor is created.
For parallel computations, an additional variable "process rank" is written into the file.
The process rank allows the user to inspect the subdomains after the computation.
The parallel `*.pvtu``*.pvtp` files are combined in a single `*.pvd` file
like in sequential simulations that can be opened with e.g. ParaView.

## Shared-memory parallelism and multi-threaded applications

Some parts of Dumux application can exploit parallelism with the shared memory model. This is for example used in the Dumux::FVAssembler by default to assemble the residual and stiffness matrix in parallel. Therefore, the assembly will usually be significantly faster on multi-processor machines.

Multithreading is enabled if a multi-threading backend is found. Currently, we support one of `OpenMP`, `TBB`,  C++ parallel algorithms, `Kokkos`. The backend is selected by `CMake` during configure and stored in the variable `DUMUX_MULTITHREADING_BACKEND`. In the `CMake` terminal output during `dunecontrol` (see @ref installation), you may see

    -- Dumux multithreading backed: TBB

or

    -- Dumux multithreading backed: Serial

if no suitable backend could be found. You can switch backends after DuMu<sup>x</sup> has been configured by running in the build folder (e.g. `dumux/build-cmake`)

    cmake -DDUMUX_MULTITHREADING_BACKEND=OpenMP .

CMake will throw an error in case the chosen backend is not available on your system.

### Multi-threaded assembly

If a multi-threading backend is found, the Dumux::FVAssembler runs in parallel per default. The assembly algorithm will first compute an element coloring and then assembly each color in parallel. The coloring prevents data races when writing into caches, the matrix, or the residual vectors. You will see output such as

    Colored 100 elements with 7 colors in 3.4625e-05 seconds.

You may disable multithreaded assembly via the command line
or the parameter file, e.g.

    ./test_executable -Assembly.Multithreading false


### Restricting the number of threads

Per default, the backend decides on the default number of threads
used. Usually, this will be the number of available threads on your system.
Often (and very important for running on clusters) you may want to restrict the number of threads used by a DuMu<sup>x</sup> application. You can do
this with the environment variable `DUMUX_NUM_THREADS`, e.g.

    DUMUX_NUM_THREADS=2 ./test_executable

Try running your application with different number of threads to see
how much your simulation time reduces.


## Hybrid parallel runs

You can combine shared memory and distributed memory parallelism.
However, note that when using distributed memory parallelism you
may not profit from additionally using multithreading and it might
be more efficient to use more MPI processes if more cores are available.

## Which components support parallel computing?

Not all components in DuMu<sup>x</sup> support parallel computing.
Currently, all supported grid implementations read the grid on the
main process and distribute after which degrades parallel efficiency.
Solver components do not support shared memory parallelism per default.
However, there is for example a parallel matrix operator (Dumux::ParallelMultiTypeMatrixAdapter) and parallel smoothers (Dumux::ParMTSSOR, Dumux::ParMTJac) available that can be used to build custom solvers.
The multidomain framework and simulations using the multidomain framework do not generally support parallel computing at the moment.
Some coupling managers (e.g. Dumux::Embedded1d3dCouplingManager) support multi-threaded assembly.

Some operations may require manual communication. See, for example,
the [diffusion example](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/examples/diffusion) for an example of manual communication to create a random initial solution that is consistent across all processes.
