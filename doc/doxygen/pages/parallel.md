# Running in parallel

In this section, we discuss how to run DuMux applications in parallel.
There is two types of parallelism implemented in DuMux which can be combined: shared memory parallelism (multi-threading)
and distributed memory parallelism.

[TOC]

## Distributed memory parallelism with MPI

The parallelization in DuMux is based on the model supported by DUNE which is based on Message Passing Interface (MPI) (distributed-
memory approach). The main idea behind the MPI parallelization is the concept of domain decomposition. For parallel
simulations, the computational domain is split into subdomains and one process (rank) is used to solve the local problem of each subdomain.
During the global solution process, some data exchange between the ranks/subdomains is needed.
MPI is used to send data to other ranks and to receive data from other ranks
The domain decomposition in Dune is handled by the grid managers.
The grid is partitioned and distributed on several nodes.
Most grid managers contain own domain decomposition methods to split the computational domain into subdomains.
Some grid managers also support external tools like METIS, ParMETIS, PTScotch or ZOLTAN for partitioning.

On the other hand, linear algebra types such as matrices and vectors do not know that they are in a parallel environment.
Communication is then handled by the components of the parallel solvers.
Currently, most solvers in DuMux are also usable as parallel solvers.
First, the Dumux::AMGBiCGSTABBackend, a parallel AMG-preconditioned BiCGSTAB solver.
Second, the Dumux::IstlSolverFactoryBackend, which provides a selections of different parallel solvers and preconditioners.
This backend makes it also possible to choose solver and preconditioner during runtime,
but this flexibility is achieved by the cost of an increased compile time.
In order for DuMux simulation to run in parallel, an MPI library (e.g. OpenMPI, MPICH or IntelMPI)
implementation must be installed on the system. However, not all parts of DuMux can be used in parallel.
Furthermore, we note that the parallel AMG preconditioner of dune-istl defaults
to an iterative SSOR coarse grid solver if no direct solver is found on your system.
Unfortunately, the iterative solver has a very high and hard-coded tolerance as a termination criterion,
which will not solve the coarse grid system with sufficient accuracy for typical problems in DuMux.
We therefore recommend to install one of the direct solver libraries supported by dune-istl.
This is either UMFPack contained in SuiteSparse, or SuperLU, see also the section on [External Libraries](#external-libraries).

## Prepare a parallel application

In order to switch to a parallel solver backend include the respective header

    #include <dumux/linear/amgbackend.hh> or
    #include <dumux/linear/istlsolverfactorybackend.hh>

Second, the linear solver must be switched to the parallel solver backend

    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>; or
    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>>;

The parallel instance of the linear solver has to be constructed with a
Dune::GridView object and a mapper, in order to construct the parallel index set needed for communication.

    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

When using the Dumux::IstlSolverFactoryBackend, solver and preconditioner have to be specified
in the file params.input by the parameters `LinearSolver.Type` and `LinearSolver.Preconditioner.Type`.
Possible solvers are bicgstabsolver or restartedgmressolver,
possible preconditioners ilu or gs.
Depending on the chosen solver and preconditioner additional parameters can be specified.
More information about the options
and the other available solvers and preconditioners can be found in `dune-istl`
(solvers.hh and preconditioners.hh).

## Run a parallel application

The starting procedure for parallel simulations depends on the chosen MPI library.
Most MPI implementations use the `mpirun` command

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

TODO
