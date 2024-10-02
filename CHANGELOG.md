Differences Between DuMu<sup>x</sup> 3.10 and DuMu<sup>x</sup> 3.9
=============================================

### General changes / structure

### Improvements and Enhancements

- __Grid I/O__: The vtu/vtp reader now allows to read unstructured grid from (ASCII) vtu/vtp files with e.g. dune-alugrid, dune-uggrid. Grid data can also be handled in parallel. Like for the GmshReader, the grid and data is read on rank 0 and then broadcasted for now.
- __Multidomain boundary__: A init function was added to coupling managers of free-flow porenetwork as well as free-flow porousmedium to allow for transient problems.

### Immediate interface changes not allowing/requiring a deprecation period:

### Deprecated properties/classes/functions/files, to be removed after 3.10:


Differences Between DuMu<sup>x</sup> 3.9 and DuMu<sup>x</sup> 3.8
=============================================

- __Requirements__: DuMux requires Dune >=2.9 and CMake >= 3.16. It was successfully tested with dune-mmesh 1.4 and OPM 2023.10.

### Improvements and Enhancements

- __Energy Balance__: Modified the energy balance implementation to correctly include the contribution of gravity.
- __Porous Medium Flow for 2pncmin__: Fixed a bug for the non-isothermal 2pncmin test to use the permeability of the current time and not the reference one as permeability can change over time.
- __TPFA Dispersion__: Fixed the transmissibility calculation for TPFA dispersion fluxes.
- __Darcy-Brinkman Freeflow__: Introduced a convenience function to add a Brinkman term to turn the Navier-Stokes model into the Darcy-Brinkman model. With this, a single domain can contain both free flow and flow through a porous medium.
The term is added as a source term in the problem using the new helper function `addBrinkmanTerm`.
The function uses a new spatial parameter interface implemented in the new `BrinkmanSpatialParams` class (`dumux/freeflow/spatialparams.hh`). The helper function can deal with isotropic and anisotropic permeabilites.
- __Face Centered Velocity Reconstrution__: For the FC Staggered discretization of free flow, a reconstrution method has been implemented to collect the full velocity vector at the face center.
- __Pore network__: Added a model and a test case for two-phase compositional fluid flow.
- __Components__: Added the new class `ShomateMethod` for calculating the heat capacity and enthalpy of components at specified temperatures. An example implementation can be found for methane (CH4).
- __SimpleH2O__: Fixed an issue where the function vaporizationEnthalpy returned the result in the wrong unit, added new nonisothermal test.
- __Facet-Coupling__:
    * Fixed the handling of duplicate degrees of freedom in the box facet-coupling model in the corner case that an internal fracture turns into a boundary fracture (see [merge request 3748](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/merge_requests/3748) for images).
    * The facet-coupling framework has been modified such that lower-dimensional domains coinciding with the bulk domain boundary are supported.
- __ISTL solver factory__: Fixed an issue with exceptions where an exception thrown during solver construction would lead to a deadlock in parallel simulation. The solver factory now communicated the failure which makes the exception recoverable, e.g. in the Newton solver.
- __Periodic Boundaries__: Fixed an issue for vector-valued unknowns. Other schemes that provide a periodic map at boundaries now also support periodicity.
- __TimeLoop__: Added function to allow adjustments of the time step size in the main file, fixes corner cases.
- __Pipeline__: Variable docker image testing is now possible directly out of dumux. A weekly Ubuntu 24.04 pipeline is introduced for a smooth transition to 24.04 in the future.
- __Dependencies__: Fieldcompare and Pylint have been upgraded to newer versions.

### Immediate interface changes not allowing/requiring a deprecation period:

- __RichardsNewtonSolver__: It is now possible to select the `MPICommunicator` used by the `RichardsNewtonSolver` (either real or dummy communicator).
- __CompositionalFluidState__: `setRelativeHumidity` was removed. Use the other setters. This setter was removed because it was very specific with a lot of specific prerequisites not fitting the general concept of the class. It was also outdated and not used in any example or test and didn't fit the index convention used in the fluid systems anymore.
- __LocalResidual of porousmediumflow/nonisothermal and porousmediumflow/nonequilibrium/thermal__: `fluidPhaseStorage` now requires the problem. The old interface is deprecated but a corresponding function call doesn't calculate anything.

### Deprecated properties/classes/functions/files, to be removed after 3.9:

- __Periodic Map__: `periodicVertexMap` has been deprecated. Use `periodicDofMap` instead.

Differences Between DuMu<sup>x</sup> 3.8 and DuMu<sup>x</sup> 3.7
=============================================

- __Requirements__: DuMux requires Dune >=2.9 and CMake >= 3.16. It was successfully tested with dune-mmesh 1.4 and OPM 2023.10.

### General changes / structure

- __Handbook__: The PDF handbook is discontinued. The content has been merged into the online documentation

### Improvements and Enhancements

- __Pore network__: Added the dual network model proposed in Koch et al (2021) https://doi.org/10.1007/s11242-021-01602-5
- __Pore network__: Added a model for heat conduction in a solid grain network
- __Nonisothermal__: Added a specialized local residual for incompressible flows where the pressure work term cancels with the pressure-dependent part of the enthalpy so that only the internal energy remains in the advective heat flux.
- __Tpfa__: The scv/scvf interfaces no longer store the geometry (can be obtained via local view of grid geometry). This means the memory footprint and also construction performance of the scv/scvfs is greatly improved, especially in higher dimensions. (Similar changes for the CVFE methods are already implemented since the previous release.)
- __Input file__: A string expression parser was added (implementation based on the ExprTK library which is shipped with DuMux). Using the new class `FunctionFromStringExpression<numArgs, ArgType>` functions based on string expression can parsed and evaluated. The most prominent use case is parsing functions
from the parameter input file. The constructor takes the expression and an array of variable names to be replaced by function arguments. The function can then
be evaluated with the function values provided in the same order as the names where specified. An arbitrary number of arguments is supported. ExprTK support very complex expressions, see https://www.partow.net/programming/exprtk/.
- __Time loop__: Fixed a bug when the time is chosen to be negative and/or the end time is set to zero. (Before we assume time is always positive and endtime never zero.). Moreover, time loops can now be constructed from and time step sizes set with `std::chrono::duration`. Finally, corner cases that may occur when using both periodic and manual check points in the `CheckPointTimeLoop` are now properly handled.
- __Chrono__: Added the conversion function `Dumux::Chrono::toSeconds` to parse strings representing numbers with trailing unit into a `std::chrono::duration`. This allows convenient construction of durations, for instance: `const auto tEnd = Dumux::Chrono::toSeconds(getParam("TimeLoop.TEnd"));`, where `TEnd = 12ms` would be valid content in the input file.
- __Hyperelastic__: Added a hyperelastic model (large deformations) and a test (in geomechanics)
- __Property system__: Added a function `inheritsFrom` that allows to query whether a type tag is present in the inheritance hierarchy of another type tag
- __PDESolver__: The LinearPDESolver/NewtonSolver now has an `apply` method that returns a bool (true if converged) instead of throwing an exception
- __Freeflow Nonisothermal__: An issue with the compositional heat flux's upwinding calculations has been fixed
- __Compositional Staggered__: The staggered discretization method has been updated to fix inconsistencies in handling density when evaluating diffusive boundary fluxes. In addition, a total mass balance is now always considered, rather than a total mole balance.
- __Parameters__: The template argument for `getParam` and `getParamFromGroup` now defaults to `std::string`, such that for string parameters one can simply write `const auto param = getParam("MyParam");`
- __Compositional Freeflow FCSFV__: Compositional transport has been implemented for the overhauled face-centered staggered finite volume scheme (FCSFV) discretization method. Ported tests are available in the testing suite.
- __GridGeometry__: No longer store the corners/geometry for scv/scvfs also for cc-mpfa and staggered discretizations. This improves the memory footprint a lot. Geometries can be obtained via the local view of the grid geometry.
- __Python bindings__: Python bindings are now disabled by default if DuMux is built without custom configuration (building without any *.opts file / passing variables to CMake). The default options in `cmake.opts` still enable the Python bindings (like before release 3.8) and thus ensure that all libraries are built when `cmake.opts` is passed to dunecontrol. (This applies also to upstream modules if cmake.opts` is used to configure them.) If you follow the default installation instructions everything should stay the same. This update tries to avoid a problem that occurs when mixing static and shared libraries and tries to make it harder for users to get the wrong setup.

### Immediate interface changes not allowing/requiring a deprecation period:

- __Newton__: The Newton solver no longer supports linear solvers without a `norm` interface when computing the resisual norm is required. The linear solvers available in Dumux all have such an interface.
- __MultiDomain__: The MDTraits are now required to state the LocalResidual type in the Subdomain policy (see multidomain/traits.hh). This only affects users that created their own MDTraits and don't use the default MDTraits.
- __PDESolver__: The PDESolver interface now requires an `apply` method that returns a `bool` instead of throwing when not converged
- DuMux-specific preprocessor macros (defined in `config.h`) are now prefixed by `DUMUX_` in order to avoid name conflicts when being passed down to other modules. Corresponding CMake variables (if existing) are now also prefixed. Explicitly, the following macros have been renamed:
  - `HAVE_KOKKOS` to `DUMUX_HAVE_KOKKOS`
  - `HAVE_OPENMP` to `DUMUX_HAVE_OPENMP`
  - `HAVE_GNUPLOT` to `DUMUX_HAVE_GNUPLOT`
  - `HAVE_CPP_PARALLEL_ALGORITHMS` to `DUMUX_HAVE_CPP_PARALLEL_ALGORITHMS`
  - `HAVE_GMSH` to `DUMUX_HAVE_GMSH`
  - `HAVE_GSTAT` to `DUMUX_HAVE_GSTAT`
  - `HAVE_PVPYTHON` to `DUMUX_HAVE_PVPYTHON`

### Deprecated properties/classes/functions/files, to be removed after 3.8:
- __MPFA__: `scvf.corners()/facetCorner()/vertexCorner()` have been deprecated. Use functions from the local view of the grid geometry instead.


Differences Between DuMu<sup>x</sup> 3.7 and DuMu<sup>x</sup> 3.6
=============================================

- __Requirements__: DuMux requires Dune >=2.9 and CMake >= 3.14. It was successfully tested with dune-mmesh 1.4 and OPM 2022.10.

### General changes / structure

- __Doxygen__: The theme of the Doxygen documentation page has been updated to a more modern look and the content has been restructured. The installation guide, the property system description and the chapter on examples and tutorials have been moved from the handbook to the Doxygen documentation. A chapter on parallelism has been added. One of the CI pipelines now includes a build of the Doxygen documentation. The result of this build can be downloaded from the job artifacts.

- __License__: DuMux is now REUSE compliant. Many files are individually licensed, others covered by rules in `.reuse/dep5`.

- __Testing/CI__: One of the CI pipelines now runs the static code analyzer `cppcheck` on every merge request.

- __Testing/CI__: The example documentation is now re-generated in the CI in order to verify that the generated README files are in sync with the sources from which they are produced.

- __Testing/CI__: The test suite now uses the [`fieldcompare`](https://pypi.org/project/fieldcompare/) library to compare VTK and data files in regression tests. The legacy comparison backend is kept as fallback.

- __Metadata__: Added a codemeta.json file describing DuMux

- __CMake__: DuMux now requires at least CMake version 3.14


### Improvements and Enhancements

- __Linear solvers__: New linear solvers have been added in `istlsolvers.hh`. They offer additional runtime options and most of them can be utilized in MPI parallel computation, eliminating the need to use `istlsolverfactory.hh` which takes longer to compile.

- __Linear solvers__: Added three multi-threaded smoothers that can be used to speedup AMG (ParMTJac, ParMTSOR, ParMTSSOR).

- __Linear solvers__: Added an iterative solver for the Stokes problem with a block-diagonal or block-triagonal preconditioner.

- __Examples__: There are now three more examples. The first example shows how to simulate diffusion using a custom model equation, with finite volume/element methods. The second one models the separation of two phases using Cahn-Hilliard model with two governing nonlinear equations. The third example models tracer spread in blood and tissue using a multi-domain model coupling 1D advection-diffusion equation and 3D diffusion equation in the porous medium.

- __Box/Diamond/Bubble/CVFE Assembly__: All CVFE schemes now use the same element solution and
the same assembler. The old assemblers have been deleted (see below).

- __IO/RasterImageWriter__: A tool now exists for writing `.pbm` and `.pgm` image files.

- __Shallow water equations__: Added new friction law `FrictionLawViscousNoSlip` for viscous thin film flow.

- __Shallow water equations__: Make regularization of the water height / roughness height optional in the friction laws.

- __Projection__: In addition to the L2-projector projecting between different grids added a helper that computes the L2-projection of analytic functions in to discrete FEM spaces (requires `dune-functions`).

- __Box__: No longer store the corners/geometry for scv/scvfs. This improves the memory footprint a lot. Geometries can be obtained via the local view of the grid geometry. (Implementation for other discretization schemes is planned for the coming release and scv/scvfs interfaces are deprecated.)

- __GridGeometry__: The grid geometry inherits from a basic discretization-agnostic implementation that can be swapped out via traits. The basic implementation contains the bounding box tree and mappers and the instance can be shared among multiple grid geometry instances reducing the number of trees/mappers that have to be instantiated for multidomain problems where both domains use the same grid, e.g. Stokes.

- __Fmt__: shipped fmt has been updated to 9.1.0. If the standard library support <format> we now use the standard library instead of fmt. Note that the standard library implementation does not support all features of fmt. In case you have been using such special features you might need to face errors for this reason.

- __Assembler__: `FVAssembler::assembleJacobian` was fixed (didn't actually compile before) and is also tested now. It assembles the Jacobian matrix. Most commonly the method `FVAssembler::assembleJacobianAndResidual` is used which also assembles the residual which is need when computing finite difference approximations of the Jacobian anyway.

- __Poromechanics__: Fixed a bug in `PoroElasticLocalResidual`, where the average density between fluid and solid was computed incorrectly, potentially leading to unphysical body forces.

- __Box/CVFE/Porenetwork Assembly bugfix__: The flux variables caches are now updated when computing the Jacobian. Before this fix, fluxvariable caches depending on the solution were not updates for CVFE schemes such that the resulting Jacobian was only an approximation.

- __Time loop__: Only insert duplicate check points once.

### Immediate interface changes not allowing/requiring a deprecation period:

- __Assembler/Newton/PDE/Solver__: We now distinguish between `SolutionVector` and `ResidualType`/`ResidualVector`. The former
contain `Dumux`-specific types like `SwitchablePrimaryVariables` as block types. The latter is native to the linear algebra backend.
The only supported linear algebra backend at the moment is Dune (common/istl). If you have implemented your own assembler or
PDE solver, you may also now need to follow this distinction in the assembler and solver interfaces. Moreover, you may need to
use specialized assign and numeric operations in case your code allows for custom block types. If you are using the classes
from the `Dumux` namespace, no change should be necessary in user code in the majority of cases.

- __Box/CVFE/Porenetwork FluxVariablesCache__: The `FluxVariablesCache` class used with CVFE method are now required to have an interface
variable `isSolDependent` that states whether the cache depends on the solution. This is used in the assembler to correctly update the cache. The `ElementFluxVariablesCache` classes with caching disabled now need a method `update` that updates the caches on demand. The
`GridFluxVariablesCache` classes with caching enabled now need a method `updateElement` equivalent to `update` for the corresponding localView for caching disabled.

- __Shallow water friction laws__: The friction laws after Manning and Nikuradse do no longer apply a limiter for the water depth. Previously, based on some application-specific assumption a roughness height was calculated that was added to the water depth. This
height can now be provided as a argument to the respective friction law class constructors
but defaults to `0.0`.

- __NewtonConvergenceWriter__: The convergence writer now takes three template arguments.
The new and last argument is the `ResidualType` (see above).

- __FaceCenteredDiamond__: The boundary treatment now follows the other control-volume
finite element schemes. This means that problem's `boundaryTypes` and `dirichlet` interfaces
are called with a sub-control volume as argument instead of
a sub-control volume face. The integration point of the scvf corresponds
to the dofPosition of the scv. In case you have been using the
`FaceCenteredDiamond` discretization and not the `boundaryTypesAtPos` and `dirichletAtPos`
version you will face an exception. Replace the functions by the version with an `scv`.
The change is made in the attempt to unify assembly for CVFE schemes

- __Box/Diamond/Bubble/CVFE Assembly__: The classes `BoxLocalAssembler`, `FaceCenteredDiamondLocalAssembler`, `PQ1BubbleLocalAssembler`, `BoxSubdomainLocalAssembler`, `FaceCenteredDiamondSubdomainLocalAssembler`, `PQ1BubbleSubdomainLocalAssembler` and corresponding
headers have been replaced by `CVFELocalAssembler` and `CVFESubdomainLocalAssembler`.

- __Default CO2 Table__: The new header `defaultco2table.hh` provides the CO2 density and enthalpy values with in the range of 290 K to 340 K and 100 KPa to 100 MPa. The script `make_co2_table.py` now generates the `co2table.hh` for a specified range, which can be directly included in the properties file.

### Deprecated properties/classes/functions/files, to be removed after 3.7:

- __AMGBackend__: `AMGBiCGSTABBackend` have been deprecated, use `AMGBiCGSTABIstlSolver` instead
- __IstlSolverFactoryBackend__: `IstlSolverFactoryBackend` now require an additional parameter `LinearAlgebraTraits`
- __BrineCO2__: Fluidsystem `BrineCO2` and binary coefficients `Brine_CO2` receiving a `CO2Table` as
template parameter has been deprecated. Pass a component instead, use `Components::CO2<CO2Table>` to
keep using a tabulated component.
- __GridGeometry__: scv/f.corner/geometry() interfaces


Differences Between DuMu<sup>x</sup> 3.6 and DuMu<sup>x</sup> 3.5
=============================================

### General changes / structure

- __Testing/CI__: Dumux is now tested with C++20 flag enabled in addition to C++17 (which is still the minimum requirement)

- __Testing/CI__: The CI system now checks for common spelling mistakes using the `codespell` tool.

### Improvements and Enhancements

- __bin/extract_as_new_module.py__: Allow for main branch named `main` or `master`

- __Components__: Fixed a bug in `TabularizedComponent` that caused data races in multithreaded applications

- __Discretization__: There is now defaults for `FluxVariablesCache` and `FluxVariablesCacheFiller` for box/cctpfa/ccmpfa/staggered
that implement an empty cache (i.e. nothing is cached for the flux variables).

- __Discretization__: Added a new face-centered FV scheme based on non-conforming FE spaces `FCDiamond`

- __Discretization__: Added a new FV scheme based on conforming FE spaces `pq1bubble`

- __Discretization__: Added a `basicGridGeometry` which isn't discretization specific and may be shared among multiple `gridGeometry`s on the same grid.

- __Grid__: Added stamped subgrids which allow for domains to be generated from repeated images

- __Parallelization__: Update grid caches in parallel

- __Parallelization__: GridView confirms whether multithreaded iteration of the grid is permitted

- __Multidomain__: Added a parallel scalar product for multidomain problems

- __Assembly__: Enabled multithreaded assembly for embedded problems

- __Multiphase and multicomponent__: Added an initial helper to choose an appropriate constraintsolver based on present phases

- __Properties__: There is now a `GetPropOr` helper that evaluates to the property type if that type is specialized for the given TypeTag and a given type if not.

- __1d3d__: Fixed a bug in the extended source stencil which didn't respect a user-defined epsilon for numeric differentiation

- __Pore-network model__: Pore-network model will no longer prevent non-wetting fluid flowing out by default. Throats  blocking non-wetting fluid can be specified by setting runtime parameter `InvasionState.BlockNonwettingPhaseAtThroatLabel`.

- __Pore network model__: Fixed a bug in the calculation of pc snap-off. To calculate the pc snap-off, corner half angle of the throat cross section shape is needed. Previous implementation was only based on square cross section. We included corner half angle in the calculation to be able to calculate pc snap-off for other cross sections than square.

- __Examples__: Extend the porenetwork_upscaling example to include non-creeping flow simulation in pore network. The example is able now to provide not only upscaled Darcy permeability but also Forchheimer permeability and coefficient (employed in Forchheimer equation).

- __Json__: Dumux now ships a fixed version of `https://github.com/nlohmann/json` to read and write JSON files. Only the symbols in the namespace `Dumux::Json` are meant to be used with Dumux. In particular, this now exports the json tree structure. More features may be added in the future.

### Immediate interface changes not allowing/requiring a deprecation period:

- __Experimental features__: All headers with experimental features have been moved into the folder `dumux/experimental`.
  Correspondingly, all unit tests for those features have been moved to `test/experimental`.

- __Stokes/Darcy__: In the coupling manager, the diffusion coefficient coupling type had a mode `Arithmetic` misspelled
  with an additional letter. The spelling has been corrected. The incorrect spelling will lead now lead to an error.

- __Box__: The grid geometry now passes information to the local view via an internal cache (instead of passing a pointer to itself directly). This allows to hide some internal interfaces. This only concerns implementers of FVElementGeometry classes that inherit from the box implementation and don't overload the constructor.
Then code possibly fails to compile. The fix is to implement the same caching concept. The cache can be a simple wrapper around the grid geometry pointer.

### Deprecated properties/classes/functions/files, to be removed after 3.6:

- __Box__: `scv/scvf.geometry()` have been deprecated, use `fvGeometry.geometry(scv/scvf)` instead
- __Box__: `scv/scvf.corner(i)` have been deprecated, use `fvGeometry.geometry(scv/scvf).corner(i)` instead
- __Cell centered__: The `computeTpfa/MpfaTransmissibilities()` interfaces now require an additional parameter `fvGeometry`
- __Staggered__: `fluxVars.advectivFluxForCellCenter()/inflowOutflowBoundaryFlux()` interfaces now require parameter `fvGeometry` instead of `element`
- __Discretization__: `Extrusion::area/volume(scvf/scv)` have been deprecated, use `Extrusion::area/volume(fvGeometry, scvf/scv)` instead
- __Richards__: Using the extended Richards model accounting for vapor diffusion in the gaseous phase by setting the property `EnableWaterDiffusionInAir` and the use of the corresponding template parameter in `richards/model.hh` has been deprecated. Use the new model `ExtendedRichards` instead
- __Navier-Stokes__: The unified `NavierStokesParentProblem` covering both momentum and mass problems has been deprecated, use separated problems instead
- __Box__: `BoxElementVolumeVariables` and `BoxGridVolumeVariables` have been deprecated, use unified `CVFE` volume variables instead

### New experimental features (possibly subject to backwards-incompatible changes in the future)

- __Meta data__: There is support for extraction of meta data in `dumux/common/metadata.hh`. The extracted features, might change in the future
and we can't give any guarantee that the name of keys and values in the extracted meta data tree will be stable between releases.

Differences Between DuMu<sup>x</sup> 3.5 and DuMu<sup>x</sup> 3.4
=============================================
- __Requirements__: DuMux requires Dune >=2.8 and CMake >= 3.13. It was successfully tested with OPM 2021.10 (which in turn requires Dune <= 2.8), see also [patches](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/patches).

### General changes / structure

- __Sequential__: The old "sequential", "IMPES-style" models have not received updates in several years,
  are separate from the rest of the code base and increasingly became a maintenance burden. Deprecation warnings
  are emitted when using the sequential models since the DuMu<sup>x</sup> 3.0 release. For DuMu<sup>x</sup> 3.5,
  these models are now removed from the code base. To continue using them use `releases/3.4`. IMPES-style algorithms
  are possible in the DuMu<sup>x</sup> 3.0 style of implementing models but their development did
  not receive much attention. See https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/merge_requests/2629
  for the state of current efforts. Further improvements and contributions (MRs) on this front
  are highly welcome (get in contact via the DuMu<sup>x</sup> mailing list).

- We adopted a code of conduct, see `dumux/CODE_OF_CONDUCT.md`

### Improvements and Enhancements

- __ParallelFor and multithreaded assembly__: For single domain applications using Box of CCTpfa discretization,
  the possibility for multithreaded (thread-parallel) assembly has been added. It is using the newly added
  `Dumux::parallelFor` interface. This features requires a backend to be available.
  This can be either external backends (TBB, Kokkos) or a working C++ STL parallel algorithms setup.
  The backend will be selected automatically, if found. You can also specify the backend by setting the
  compiler definition `DUMUX_MULTITHREADING_BACKEND=TBB,Cpp,Kokkos,Serial`, where Serial forces to always run in serial.
  For the assembly, you can explicitly turn multithreading off setting the runtime parameter `Assembly.Multithreading = false`.
  If a backend is available and the discretization allows it, the default is multithreaded assembly.
  When using TBB or Kokkos is it required (recommended) to use `Dumux::initialize` in the main file.

- __initialize__: A new control function `Dumux::initialize` has been added which initializes shared and distributed
  memory parallelism helpers. It is recommended (and may be required for multithreaded applications) to use `Dumux::initialize`
  instead of manually initializing the `Dune::MPIHelper`.

- __Problem/Spatialparams interface__: The interface functions `extrusionFactor`/`extrusionFactorAtPos` and `temperature`/`temperatureAtPos` have been removed from the problem interface and added to the spatial params interface. There is now a default set
for temperature (`293.15K`) which can be changed via the input file parameter `SpatialParams.Temperature`. For more
flexibility overload the function in your spatial parameter class like usual.

- __material/spatialparams__: The folder has been dissolved and the headers have been moved to other folders. For example,
spatial parameter classes specific to porous medium flow problems have been moved to `dumux/porousmediumflow`. This is due
to the realization that spatial parameters do more than defining the "material" (which refers to the porous medium).
For example, free flow problem also have parameters that may spatially vary (such as the extrusion factor in an extruded domain).

- __Tensor average__: Added a function `faceTensorAverage` which performs a particular average of a tensor at interfaces.
This replaces the interface `harmonicMean` in the spatial parameters.

- __Component__: `gasViscosityIsConstant` added to component interface

- __Fluidsystems__: Implemented `viscosityIsConstant(int phaseIdx)` for `H2OAir`

- __Dispersion__: Dispersion fluxes have been added as an option for the compositional and thermal one-phase porous medium flow models. These models use either a Scheidegger-type dispersion tensor, which is dependent on the velocity field and two length parameters, or a full and constant (_not_ solution-dependent, but possibly spatially varying) tensor that can be user defined in the spatial parameters. For compositional models coupled with flow models (e.g. 1pnc), using the Scheidegger-type dispersion tensor is only implemented for the box discretization method.

    To enable either thermal or compositional dispersion, please define these properties within your `properties.hh` header. For example:
    ```cpp
    template<class TypeTag>
    struct EnableCompositionalDispersion<TypeTag, TTag::MyTest> { static constexpr bool value = true; };
    template<class TypeTag>
    struct EnableThermalDispersion<TypeTag, TTag::MyTest> { static constexpr bool value = true; };
    ```

    To determine the type of dispersion tensor, please define the property `CompositionalDispersionModel` within your `properties.hh` header. Per default, the `ThermalDispersionModel` is set to the same model as the `CompositionalDispersionModel`, but this can be set specifically as well. For example:
    ```cpp
    template<class TypeTag>
    struct CompositionalDispersionModel<TypeTag, TTag::MyTest> { using type = ScheideggersDispersionTensor<TypeTag>; };
    template<class TypeTag>
    struct ThermalDispersionModel<TypeTag, TTag::MyTest> { using type = FullDispersionTensor<TypeTag>; };
    ```

    The parameters describing your dispersion tensor can then be included in your `spatialparameters.hh` file, and passed via input parameters. An example of this can be seen in the `test/porousmediumflow/1pnc/dispersion/` folder, and in the `test/porousmediumflow/tracer/constvel/` folders.

- __Embedded coupling__: Add a coupling manager for the 1D-3D projection based scheme with resolved interface introduced in Koch 2022 (
https://doi.org/10.1016/j.camwa.2022.01.021)

- __RANS Boundary Types__: Wall boundaries for the RANS turbulence models can now be set using the `setWall` method in the `RANSBoundaryTypes` class. This replaces the old `isOnWall` and `isOnWallAtPos` functions.

- __Discretization tags__: We introduced tags in the namespace `DiscretizationMethods` (with s) for each discretization method.
  These tags replace the `enum class DiscretizationMethod`. Tags have several advantages over the enum. Each tag is a named type
  (see `dumux/common/tag.hh`) so they can for example be used in tag dispatch. Moreover specializing with tags is extensible.
  For example we specialize the template `DarcysLawImplementation` for each discretization. When using the enum no new discretization
  methods can be added without changing `enum class DiscretizationMethod`. With tags you can make your own tag and specialize a
  class for that tag. This means new discretization methods can be developed in a modular fashion. The introduction of tags
  involves a couple of non-backwards-compatible changes, mostly in implementation classes that shouldn't affect most users (see below).

- __Box__: The box scheme now supports the case that volume variables depend on all dofs in the element. In that case, previously
  only a Jacobian approximation was assembled. As computing the added derivatives causes a significant
  overhead in the assembly, the feature is disabled per default. To enable the full Jacobian in this case set the parameter
  `Assembly.BoxVolVarsDependOnAllElementDofs = true`. The new feature is tested in the new test `test_2p3c_surfactant` which
  features relative permeability that depend on the pressure gradient. The test only passes with the full Jacobian.

- __Local views__: The bind function associated with the local view of the FVElementGeometry, the ElementVolumeVariables, and the ElementFluxVariablesCache have been simplified. Now it is possible to directly create each of these objects where they are already bound. The following are examples of each new call:

    ```cpp
    const auto fvGeometry = localView(gridGeometry).bind(element);
    const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);
    const auto elemFluxVarsCache = localView(gridFluxVarsCache).bind(element, fvGeometry, elemVolVars);
    ```
    This is also available for the `bind()` `bindElement()` and `bindScvf()` functions. The existing methods for binding will remain.

    Please note however that when working with element loops, separating construction and binding of the local view is more efficient, i.e.

    ```cpp
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
        fvGeometry.bind(element);
    ```

- __Construction and update of GridGeometries changed__: Grid geometries are fully updated after construction. Additional call of update functions are therefore only needed after grid adaption. Calling the update functions after construction now leads to a performance penalty.

- __Geometry__:
    - Added implementation of sphere and bounding sphere approximation algorithms
    - Added distance queries for Point->BoundingBoxTree
    - Added DistanceField - a wrapper around a geometry set for fast distance queries
    - Added WallDistance - a utility to compute distances to the domain boundary (e.g. for RANS wall functions)
    - Added 3D-3D intersections
    - Added 2D-2D intersections in 3D
    - Fixed a wrong epsilon for floating comparisons in the 2D-2D intersection algorithm

- __Parallel grid partitioning__: Added a simple backend for using Scotch as partitioner for a grid read on one process (e.g. when using the Dune Gmsh or DGF readers). Repartitioning is not implemented yet.

- __Cake grid creator__: The cake grid creator can now be used in parallel simulations

- __Privarswitch__: Fixed a bug in the privar switch which did not fully reset the `switched` variable. This lead
to a possibly increased number of Newton iterations.

- __1pnc__: `delp` has been removed from the default vtk output fields in an attempt to streamline output between models.
All information present in `delp` is present in `p` and `delp` has generally a reduced precision that also leads to some
problems in automated testing.

- __Richards__: the `Richards` model now works together with the generic `FluidSystem::TwoPImmiscible` as long as a gas
and a liquid phase are present.

- __Richards__: Fixed a bug that creeped into the 1.5-phase model so it actually computes diffusion in the gas phase now.

- __Richards__: Fixed a bug that checked that the fluid viscosity is _not_ constant whereas the check should have asserted that it is constant

- Fixed `intersect` in `SplineCommon`

### Immediate interface changes not allowing/requiring a deprecation period:

- __Discretization tags__: The following classes have changed from a template argument of type `DiscretizationMethod` (an `enum`) to
  class type template arguments and are now specialized for tags: `TwoPScvSaturationReconstruction`, `ScvfToScvBoundaryTypes`, `ProblemTraits`, `FluxStencil`, `EffectiveStressLaw`, `HookesLaw`, `FacetCouplingManager`, `FacetCouplingMapper`. Moreover this affects the following helper classes: `IsFicksLaw`, `CheckOverlapSize`, `PartialReassemblerEngine`, `SubDomainAssemblerType`. Finally, this change has been made for many implementation classes. These should not have been used directly anyway, so we do not explicitly list them here. Examples are `DarcysLawImplementation`, `LinearSolverTraitsImpl`. See !2844 for more details.
  If you face a compiler error from these changes, contact us. Most like the solution is to simply try replacing occurrences of `DiscretizationMethod` with the corresponding tag from `DiscretizationMethods`, or types `DiscretizationMethod` in template arguments by generic `class`.

- __Coupling managers__: The coupling managers now store shared_ptrs of the subdomain solutions to be able to manage memory outside. There is compatibility interface that is deprecated but it won't allow for assignments
  of the form `this->curSol() = sol;` since curSol returns a MultiTypeVector of references. If assignment is needed for some reason, use a hybrid loop over all subdomain indices.

- __Virtual interface of GridDataTransfer__: The `GridDataTransfer` abstract base class now required the Grid type as a template argument. Furthermore, the `store` and `reconstruct` interface functions do now expect the grid as a function argument. This allows to correctly update grid geometries and corresponding mapper (see "Construction and update of GridGeometries changed" above in the changelog)

- `PengRobinsonMixture::computeMolarVolumes` has been removed without deprecation. It was used nowhere and did not translate.

- __ShallowWaterViscousFlux__: The template argument `PrimaryVariables` has been removed without deprecation. It was not needed.

### Deprecated properties/classes/functions/files, to be removed after 3.5:

- `update()` functions of grid geometries, which do not receive the `gridView`, are deprecated, use `update(gridView)` instead.
- `enum class DiscretizationMethod` and associated functions, to be replaced by tags
- `test_dumux.sh` is deprecated.
- `compareparameters.sh` is deprecated, use `generate_parameterlist.py` instead.
- `replace_property_macros.sh` is removed.
- `isOnWallAtPos` and `isOnWall` are no longer used in the RANS models to designate wall boundaries. These boundaries are now directly set in the RANSBoundaryTypes using the setWall() function.
- the `temperature` and `extrusionFactor` interfaces in the problem class have been deprecated and have been moved to the spatial parameters.
- Porous medium flow models should now inherit from the new base spatial parameters that can be found in the folder `dumux/porousmediumflow/`, which allow users to overload the new `temperature` and `extrusionFactor` interfaces.
- Free flow and pore network models now also expect the user problems to expose spatial parameters, in which `gravity`, `temperature` and `extrusionFactor` are defined. The respective problem interfaces have been deprecated.
- `harmonicMean` has been deprecated in the spatial params use new `faceTensorAverage`
- The problem interfaces for fluid properties in the poroelastic model, namely `effectiveFluidDensity` and `effectivePorePressure`, have been deprecated and were moved to the spatial parameters.
- The function `shearStress` in the class `FrictionLaw` and its child classes has been renamed to `bottomShearStress` and the return value has been changed. The function returns now the actual bottom shear stress and not the bottom shear stress term as it used in the shallow water model. The bottom shear stress term of the shallow water model is the bottom shear stress multiplied with minus one and normalised by the water density.

### New experimental features (possibly subject to backwards-incompatible changes in the future)

- __Staggered grid__: The staggered grid implementation has been overhauled. Unfortunately, this overhaul has not been completed yet.
Most of the Navier-Stokes tests now use the new implementation. The old implementation is still available and not deprecated yet,
but will be phased out after the release. For now both implementation live next to each other in the code base.
The new implementation is more closely built on the multidomain framework and now fully
realizes the finite volume abstractions described in the Dumux paper. The first aspect means that mass and momentum are now
separate sub-models of Navier-Stokes than can be (and are) discretized with different discretization methods and then coupled
together via coupling managers. The implemented mass discretization is CCTpfa. The momentum discretization is a face-centered
staggered finite volume scheme (FCSFV). The second aspect means that for the FCSFV scheme, subcontrol volumes and faces
are now represented by corresponding classes in the code just like for CCTpfa, Box, CCMpfa.
There is a problem class added that helps to implement both the mass and the momentum problem in one (templated) class.
Boundary conditions are now clearly separated into mass and momentum boundary conditions.
When the new implementation is fully adapted the documentation will be updated with it, this might take some time
and is not included in this release yet.

- __FF-PNM__: Added a model and test for coupling a porenetwork with a freeflow domain (see Weishaupt PhD: http://dx.doi.org/10.18419/opus-10932)

- __Python bindings__: The Python bindings work with Dune master now, which features an improved way of installing them.
The new way will be described better in the documentation once Dune 2.9 is release. Until then we refer to the
documentation of Dune. The setup with Dune 2.9 is not compatible with the setup with Dune 2.8 but we made sure
that Dumux 3.5 support both variants.

- __Pore-network models__: The development continues and many smaller things have been improved.
The PNM models remain experimental. The grid creator has been improved in usability. Added a
convenience script to extract PNM with porespy and create a DGF grid usable with Dumux.

### Continuous integration

- __Python bindings__: The Python bindings are now tested in the CI. Furthermore, Python scripts are
  automatically checked with linters, see `dumux/python/README.md` for more information.


Differences Between DuMu<sup>x</sup> 3.4 and DuMu<sup>x</sup> 3.3
=============================================

### Improvements and Enhancements
- __Requirements__: DuMux still requires Dune >=2.7 and CMake >= 3.13. It was successfully tested with OPM 2021.04.

- __Pore-network models added to DuMux__:
    - Three fully implicit pore-network models (1p, 1pnc, 2p) have been added.
    - A quasi-static 2p PNM for the creation of pc-Sw curves has been added.
    - A new namespace `Dumux::PoreNetwork` has been introduced, containing all relevant classes and functions for pore-network modeling.
    - An example introduces the use of the 1p PNM for the estimation of the upscaled Darcy permeability.
    - Advection type including inertia effect for simulations in the non-creeping flow regime is available.
    - Note that this is still considered a rather _experimental_ feature. Everything within namespace `Dumux::PoreNetwork` might undergo (backward-compatibility breaking) changes _without prior notice_.

- __Several scripts have been translated to Python__:
    - `installexternal.sh` to install external modules and libraries now rewritten as python script `bin/installexternal.py`
    - `getusedversions.sh` to extract the used DuMux/Dune versions of a module (new script: `bin/util/getusedversions.py`)
    - `extractmodulepart.sh` no longer creates an install file, instead, you can now generate install scripts for your module using the new script `bin/util/makeinstallscript.py`.
    - Note: the old shell scripts will be removed after release 3.4.

- __Python bindings__: There is now a small finite volume example included in the Python tests using the wrapped grid geometry and problem, see [test_explicit_transport_cctpfa.py](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/releases/3.4/test/python/test_explicit_transport_cctpfa.py).

- __Law for water vapor viscosity in gas mixtures changed__: Polynomial laws for determining the viscosity of vater vapor in gas mixtures are introduced and applied by the new function `h2oGasViscosityInMixture`. The polynomial laws give a better approximation of the gas viscosity especially at higher temperatures ( >273.15 K) and a high water vapor content.

- __Newton line search__: The line search algorithm decreases the step size until the residual decreases. The lower bound
of the step size can now be set in the input file via the parameter `Newton.LineSearchMinRelaxationFactor`.

- __Material / Constant component__: The `Component::Constant` can now be used in non-isothermal simulation. Simple relations
for internal energy and enthalpy depending on temperature and constant heat capacity have been added.

- __Linear PDE solver__: The `LinearPDESolver` can reuse the matrix and thus avoid unnecessary reassembly. See [test/porousmediumflow/tracer/constvel/main.cc](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/releases/3.4/test/porousmediumflow/tracer/constvel/main.cc#L119) for an example.

- __Ordering strategies for UMFPack__:
It is now possible to [choose an ordering strategy](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/releases/3.4/dumux/linear/seqsolverbackend.hh#L851) for UMFPack via the runtime parameter `LinearSolver.UMFPackOrdering` or by calling the `setOrdering()` method of `UMFPackBackend`. This can have a positive effect on the solver's performance, depending on the matrix structure.

- __Add setRetryTimeStepReductionFactor() function  to NewtonSolver__:
This function allows to set the factor by which the time step is reduced after a failed iteration. Can be used, e.g., for custom Newton solvers inheriting from this class and using a more sophisticated time step management.

- __Improve upwind interface__:
    - Additional functions have been added to `upwindscheme.hh` which can be used to `apply` the upwind scheme without the need of `FluxVariables` and to get the `multiplier` by which the flux is multiplied.
    - These additional functions are now used in `CCTpfaForchheimersLaw` by providing the upwind scheme as additional template argument.

- __Simplified use of SubGridManger__: It is now possible to specify pixel dimensions when reading an image file used as grid. For instance, `Grid.PixelDimension = 1e-3 1e-3` will scale the domain automatically such that the grid cells have a side length of 1e-3 m and you don't need to specify `Grid.UpperRight` anymore.

- __String utility__: There is a new header `common/stringutilities.hh` that implements two functions `tokenize` and `split`
that can split strings at a given delimiter.

- __Add linearRegression() function  to math.hh__:
This function gets a set of (x, y) data and calculates/returns the intercept and the slope of the regression line using the standard least squares method.

- __Shallow water__: Added a heuristic turbulence model based on a mixing length and resulting in a turbulent viscosity term.

### Immediate interface changes not allowing/requiring a deprecation period:
- __Newton__: The global parameter defaults have been substituted for local parameter defaults (in nonlinear/newtonsolver.hh). If
              you have been relying on global defaults (reading parameters without supplying a value in the input file nor a default)
              you will get a runtime ParameterException. To solve this simply provide a default or set the value in the input file.
- __MPNC__: The `MPAdapter` can now also be called with a temporary `pcKrSw` object. For this, the compiler needs to deduce the
            class' template argument types. You may need to adapt your `spatialParams` from
```
using MPAdapter = Dumux::FluidMatrix::MPAdapter<PcKrSwCurve, 2>;
...
auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
{
    return makeFluidMatrixInteraction(MPAdapter(pcKrSwCurve_));
}
```
to
```
// alias for MPAdapter is removed
auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
{
    return makeFluidMatrixInteraction(FluidMatrix::MPAdapter(pcKrSwCurve_));
}
```
- __Grid geometry__: The local views of grid geometries are now required to implement the interfaces
`element()` (returning the bound element) and `isBound()` returning a `bool` which is `true` if the
functions `bind` or `bindElement` have been called (i.e. the local geometry is in a bound state). These
interfaces are currently not used (except in the unit tests) but will be required
by the assembler in future DuMux versions.

### Deprecated properties/classes/functions/files, to be removed after 3.4:
- `Dumux::IsIndexable<T, I>` is deprecated, use `Dune::IsIndexable<T, I>`directly.

- The property `NumEqVector` has been deprecated. The class `NumEqVector` is now defined in the namespace `Dumux` in the header file `dumux/common/numeqvector.hh`.

- The member function `update()` of mappers is deprecated, use `update(gridView)` when updating the mapper after a grid or grid view change.

- All custom mapper implementations should implement `update(gridView)` replacing `update()`. Mappers with `update()` will no longer be supported after support for Dune 2.7 is dropped.

### New experimental features (possibly subject to backwards-incompatible changes in the future)

- __Time stepping__: a first implementation of a generic `Dumux::MultiStageTimeStepper` was introduced,
that allows for using different time integration schemes besides the currently supported implicit and explicit
Euler schemes. However, this poses new requirements on the linear system assemblers, which are not yet met by
the standard assemblers provided in DuMux. This will be extended in future versions.

### Continuous integration

- A first version of the DuMux GitLab CI is now at the disposal of all developers. While every night a complete
test pipeline is triggered automatically, developers have the possibility to manually start test pipelines
within merge requests that automatically identify and run only those tests affected by the changes introduced.


Differences Between DuMu<sup>x</sup> 3.3 and DuMu<sup>x</sup> 3.2
=============================================

### Improvements and Enhancements

- __Requirements__: DuMu<sup>x</sup> now requires Dune >=2.7 and CMake >= 3.13.
- __New way to use material laws__: The usage of laws for pc-Sw and kr-Sw has been completely revised. A caller does not have to pass a `parameters` object to the laws anymore. The `spatialParams` now provide a `fluidMatrixInteraction` function which bundles an arbitrary number of
 different interaction laws such as a pc-Sw and kr-Sw curve and interfacial areas.
 New pre-cached spline laws were added which can help to increase efficiency. The usage of the old interface is deprecated and warnings will be raised. The old interface will be removed after the release of 3.3.
- __New example__: We have added another free-flow example dealing with lid-driven cavity flow.
- __Install script written in Python__: The DuMu<sup>x</sup> install script has been translated to Python to improve portability. The old shell script will be removed after release 3.3.
- __Improved velocity reconstruction__: The velocity reconstruction for immiscible porous-media models has been improved, leading to slightly
  different velocity fields in the vicinity of Neumann boundaries.
- __Python bindings (experimental)__: Basic support for Python bindings has been added. Python bindings are an experimental feature
  and might undergo unannounced API changes until further notice. This concerns the files in the folders `python` and `dumux/python`. To activate
    - add `-DDUNE_ENABLE_PYTHONBINDINGS=TRUE` and `-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE` to your CMAKE_FLAGS and run dunecontrol
    - run `python3 dune-common/bin/setup-dunepy.py`
    - adapt your PYTHONPATH environment variable as described [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/python)
- __Obtain DuMux via pip__: The Dumux source code can now be installed using `pip install dumux`. This works for the Python bindings as well as the C++ source code. Matching Dune modules are also installed. Currently not all dune modules have Python packages available. We recommend trying this new feature in a virtual environment by typing
```sh
python3 -m virtualenv venv
source venv/bin/activate
pip install dumux
```
To test the setup you can use the following Python script
```py
from dune.grid import structuredGrid
from dumux.discretization import GridGeometry

gridView = structuredGrid([0,0],[1,1],[5,5])

gridGeometry = GridGeometry(gridView, discMethod="cctpfa")
gridGeometry.update()

print("The total number of scvs is {}".format(gridGeometry.numScv()))
print("The total number of scvfs is {}".format(gridGeometry.numScvf()))
```
- __fmt-library__: We now include a basic version of the [fmt-library](https://github.com/fmtlib/fmt) which implements `std::format` (coming with C++20) without the need for C++20.
  In order to use this, include `<dumux/io/format.hh>`. `format`, `format_to`, `format_to_n`, `formatted_size` are available in the `Dumux::Fmt` namespace.
  The string formatting is documented [here](https://en.cppreference.com/w/cpp/utility/format/formatter#Standard_format_specification) and follows the Python string formatting rules.
  The functions are documented on [cppreference](https://en.cppreference.com/w/cpp/utility/format).
 - __RANS__: The RANS models now include variable densities. Compositional or nonisothermal RANS models could produce slightly different, more accurate, results.

### Immediate interface changes not allowing/requiring a deprecation period:
- __Flash/Constraintsolver__: The flashes depending on material laws are immediately required to use new-style material laws (fluidMatrixInteraction interface in spatialparams)
- __Box interface solver__: The box interface solver immediately requires the new material law interface without deprecation period. Use the new class `BoxMaterialInterfaces` and update your spatial params to use the new fluidmatrixinteraction interface to be able to use the box interface solver in version 3.3.
- For the "sequential" models, the property `BoundaryTypes` has been simply renamed to `SequentialBoundaryTypes`
- __Quadmath__: Dumux::Quad has been removed without deprecation. Use Dune::Float128 instead.
- Within the RANS group, two additional runtime parameters have been included 'IsFlatWallBounded' and 'WriteFlatWallBoundedFields'.
For both the K-Epsilon and Zero-eq RANS models the 'IsFlatWallBounded' runtime parameter should be set as True,
as wall topology is not supported for these models with our geometric constraints. If not set as true, the geometry
will be checked before the model is run. If either the runtime parameter or the geometry check indicate non-flat walls,
the model will terminate. To add FlatWallBounded specific output to the vtk output, WriteFlatWallBoundedFields can be set as True.
- __1d3d coupling__: The kernel coupling manager has been replaced with the one from Koch et al (2020) JCP https://doi.org/10.1016/j.jcp.2020.109370
- __1d3d coupling__: The average and surface coupling managers has been updated with a slightly more accurate version to compute the stencils and the average operator.
The results might differ a bit when using coarse grids. However, both version are expected to converge to the same result with grid refinement.

### Deprecated properties/classes/functions/files, to be removed after 3.3:

- The property `BoundaryTypes` has been deprecated. The boundary condition type can now be deduced from the problem type using `ProblemTraits`.

### Deleted classes/files, property names, constants/enums:

- Everything that has been deprecated before release 3.2 has been removed.
- All of the geometry headers previously saved in `dumux/common/geometry` have been relocated to `dumux/geometry`.
  The headers in `dumux/common/geometry` are deprecated and will be removed after 3.3. The geometry tests have been moved from `test/common/geometry`
  and `test/common/boundingboxtree` to `test/geometry`.

### Other noteworthy changes:
- Releases earlier than 3.0 are no longer automatically tested or supported.

Differences Between DuMu<sup>x</sup> 3.2 and DuMu<sup>x</sup> 3.1
=============================================

### Improvements and Enhancements

- __C++17__: DuMu<sup>x</sup> now requires a C++ compiler supporting the C++17 features of GCC 7 (e.g. GCC 7, Clang 5).
- __Radially symmetric problems__: We now have support for radially symmetric problems (disc, ball, toroid). The support comes in form of wrappers for sub control volumes and faces that overload the respective `volume()` and `area()` function turning a 1d or 2d problem into a 2d or 3d radially symmetric problem.

- __Improvements of Beavers-Joseph(-Saffman) condition for the free flow model__: The naming for handling BJ(-S) boundary conditions has been adapted from `isBJS()` to `isBeaversJoseph()` / `setBJS()` to `setBeaversJoseph()`. In order to consider the velocity within the porous medium, the old `velocityPorousMedium(element, scvf)` method (returning a Scalar) has been renamed to `porousMediumVelocity(element, scvf)` (returning a velocity vector). The latter defaults to `VelocityVector(0.0)`.

- __Van Genuchten__: The VanGenuchten-Mualem material law now allows to set a parameter `l` (default to 0.5) which is sometimes fitted.

- __Runtime variable output precision e.g. Float64__: The VtkOutputModule has been adapted to allow easy changes of the vtk output precision. It is now possible to specify output precision in the input file using `Vtk.Precision` followed by either `Float32`, `Float64`, `UInt32`, `UInt8` or `Int32`. `Float32` stays the default. We especially advice the use of `Float64` when working with restart files.
An additional new option is `Vtk.CoordPrecision` which changes the precision of the coordinates only and uses the default of `Vtk.Precision`.

- __Effective Laws and Diffusion Coefficients__: The effective laws interface has been changed within !1684. The interface for these laws has been unified, and all coefficients are to be stored in containers that fit to the model. These quantities should then be added in the volumeVariables, meaning all effective quantities would be accessible from the volumevariables.

- __Examples__: The documentation of the examples has been improved further, focusing on readability and convenience. Further, three additional examples are included the folder `examples`. To get an overview, point your browser to https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples.

- __Linear solvers__: There is a new ISTL solver factory backend which allows to choose solvers at runtime (requires dune-istl 2.7) and also enables more parallel solvers (requires dune-istl > 2.7.0).


### Immediate interface changes not allowing/requiring a deprecation period

- Remove `Grid.HeapSize` as dune-ugrid removed the according feature as well.
- __Van Genuchten__: Corrected VanGenuchten-Mualem exponent in the nonwetting saturation formula (`1/3` instead of `1/2` (or `l`, see above))
- __Van Genuchten__: Corrected VanGenuchten-Mualem implementation of `dkrn/dSw`
- __Brooks-Corey__: Corrected Brooks-Corey implementation of `dkrn/dSw` and added the derivatives for the regularized version
- __AMGBackend__: The internal structure of the AMGBackend and the ParallelISTLHelper has been overhauled, as only used by the AMG, we did not make the changes backwards-compatible
- The global default parameters for linear solvers have been removed and moved to the class `LinearSolver`.
This only affects users that directly obtain this parameter via `getParam` somewhere in the code.

- __Sequential linear solver backends__: Remove template argument `precondBlockLevel` from `solve` functions. The preconditioner block level is now determined automatically, assuming a value of
1 for regular BCRS matrices and a value of 2 for MultiTypeBlock matrices. The respective calls from the `NewtonSolver` and `PDESolver`classes have been adapted.

- __Change matrix block arrangement for staggered models__: The matrix block structure has been adapted such to comply with the literature standard, i.e., having the velocity block (A) on `M[0][0]`
rather than on `M[1][1]`. This also requires re-arranging the submodels and properties in dumux-multidomain such that the face-related classes and vector entries now appear before the cell-centered ones.

```math
M = \begin{pmatrix}
  D & C\\
  B & A
\end{pmatrix}

\qquad => \qquad

M = \begin{pmatrix}
    A & B\\
    C & D
    \end{pmatrix}
```

Backwards-compatibility can only be provided to a certain extent. The following changes need to made in the main file:

1.) change the order of the arguments for the `assembler` such that it reads:
```c++
auto assembler = std::make_shared<Assembler>(std::make_tuple(ffProblem, ffProblem, otherProblem, ...),
                                             std::make_tuple(ffGridGeometry->faceFVGridGeometryPtr(),
                                                             ffFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                             otherFvGridGeometry, ...),
                                             std::make_tuple(ffGridVariables->faceGridVariablesPtr(),
                                                             ffGridVariables->cellCenterGridVariablesPtr(),
                                                             otherGridVariables, ...),
                                             couplingManager,
                                             timeLoop, solOld);

// Not changing the arguments will yield a deprecation warning stating this hint but the code still compiles and runs.
```

2.) change the order of arguments in the `partial` function:
```c++
ffSol = partial(sol, ffFaceIdx, ffCellCenterIdx);

// Not changing the argument will rise a compiler error which makes the MR not fully backwards-compatible.
```

Regarding changes made to the effective laws and diffusionCoefficient containers, Backwards-compatibility is maintined for a large extent, barring any volumevariable classes defined externally that inherit from the non-isothermal volumevariables.
If you use a self defined volumevariables class that inherits from the non-isothermal volumevariables, please adapt the your volumevariables class to fit to the non-isothermal volumevariables, and include the new methods for accessing the diffusion and effective diffusion coefficients.

- Tracer model: tracer fluid systems do no longer provide a `getMainComponent` function since this simply doesn't make sense -- the main bulk component is not modeled.

### Deprecated properties, to be removed after 3.2:
- __GridView__: The property `GridView` has been deprecated and can be accessed via `GridGeometry::GridView` instead.

### Deprecated classes/files, to be removed after 3.2:
- __AMGBackend__: The class AMGBackend is deprecated and has been replaced by AMGBiCGSTABBackend which gets some different template arguments
- __AMGTraits__: AMGTraits are deprecated, are to be replaced by LinearSolverTraits and restructured internally. As they were only used by the AMGBackend, we did not make the internal changes backwards-compatible

### Deprecated member functions, to be removed after 3.2:
- __DiffusionCoefficient(various arguments)__: These coefficients are now defined in the volvars and stored in a container fit to the model. To access these values, use the unified ```c++ diffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) ```
- __EffectiveDiffusivity(various arguments)__: These values are now defined in the volvars and stored in a container fit to the model. To access these values, use the unified ```c++ effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) ```

### Deleted classes/files, property names, constants/enums
Everything that has been deprecated before release 3.1 has been removed.

Differences Between DuMu<sup>x</sup> 3.1 and DuMu<sup>x</sup> 3.0
=============================================

### Improvements and Enhancements

- __Examples__: Three extensively documented examples have been added which highlight interesting features and show how to apply DuMu<sup>x</sup> to interesting problems. They can be found in the new folder `examples`. To get an overview, point your browser to https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples.

- __Porousmediumflow__: Added a new porous medium model for the energy balance of a porous solid (general heat equation).

- __Multidomain__: It is now possible to use the facet coupling module together with the mpfa-o scheme in the bulk domain.

- __Box__: The box scheme works now on grids with prism / wedge elements in 3D.

- __Diffusive fluxes__: We revised the formulation of the diffusion laws to allow for changes between mass-averaged velocity reference systems and molar-averaged velocity reference systems. The standard formulation is now set to mass-averaged velocity reference systems.

- __GridManager__:
    * Supports now reading unstructured grids and data from vtu/vtp files (ASCII, XML format) sequential.
      For UGGrid and FoamGrid you can specify a grid in such a format in the input file:
      `Grid.File = mygrid.vtu` / `Grid.File = mygrid.vtp`. The associated data is then available in the grid data object.
    * Instead of always including all grid manager specializations you can now only include the specialization that you need.
      For example, if you only use YaspGrid in your code, you only need to include `dumux/io/grid/gridmanager_yasp.hh`. The
      convenience header `dumux/io/grid/gridmanager.hh` still includes all specializations.
    * Added a `NetPBMReader` which allows to read simple raster images files (`*.pbm` and `*.pgm`). Can be used, e.g., with
      `dune-subgrid` in order to create a grid from an image file.

- __Freeflow__: A second order approximation of the convective term of all Navier-Stokes based models is now available.
  This can be enabled using the property `UpwindSchemeOrder`. This property defaults to a first order upwinding approximation,
  which is the method used for all models in release 3.0. If this property is set to `2`, a second order flux limiter method will be used.
  Various flux limiter functions have been implemented to maintain the monotonicity of this discretization. Per default the
  `Minmod` flux limiter is used, but `Vanleer`, `Vanalbada`, `Superbee`, `Umist`, `Mclimiter`, and `Wahyd` based flux limiters
  are also available. These can be specified using the input file entry `Flux.DifferencingScheme`. These methods are also
  implemented for non-uniform structured grids (e.g. YaspGrid - TensorProductCoordinates). Per default, a scheme assuming a
  uniform grid is used, but two other methods, `Li` and `Hou`, are both available for adaptations to non-uniform grids.
  These can be specified using the input file entry `Flux.TVDApproach`.

- __RANS__: The called RANS model, defined in the properties system, will specify, via the model traits,
  which RANS problem implementation should be used. In each problem file, the initial and boundary conditions can be
  set using templated functions based on the model type. Examples of these functions exist in the RANS based tests.
  No further preprocessor macros are required.

- __ShallowWater__: Thanks to Leopold Stadler we now have a shallow water equations solver / model. Have a look at freeflow/shallowwater and give it a try with the dam break test at test/freeflow/shallowwater. The are also some friction laws computing shear stresses (Manning, Nikuradse) to account for friction in a channel/river bed, thanks to Martin Utz.

- __Staggered__: Added a `StaggeredNewtonConvergenceWriter` for the staggered grid discretization scheme.

- __Solver__: There is a new abstract base class `PDESolver` that is a class that does linearize & assemble, solve and update.
  The NewtonSolver now derives from this class (interface is unchanged). A new class `LinearPDESolver` simplifies solving linear problems
  by reducing the code in the main file and streamlining the terminal output to look like the Newton output.


### Immediate interface changes not allowing/requiring a deprecation period

- `NewtonConvergenceWriter`'s first template argument has changed from `GridView` to `FVGridGeometry`. This allows to call the `resize()` method after a grid change without any arguments.
  Here is an example how to instantiate the convergence writer:
  ```
  using NewtonConvergenceWriter = Dumux::NewtonConvergenceWriter<FVGridGeometry, SolutionVector>;
  auto convergenceWriter = std::make_shared<NewtonConvergenceWriter>(*fvGridGeometry);
  ```

- The interface of the abstract `TimeLoopBase` class has been extended by `virtual void advanceTimeStep()`, `virtual void setTimeStepSize(Scalar dt)`,
  `virtual Scalar maxTimeStepSize()`, and `virtual bool finished()`, thus forcing the inheriting classes to implement those functions.
  `TimeLoop` is no longer a template argument in `NewtonSolver`'s `solve()` method, thereby increasing
  type safety and providing better compiler error messages.

- __RANS__: The base problems for all turbulence models e.g. `ZeroEqProblem` have been been renamed to one generic `RANSProblem`

### Deprecated properties, to be removed after 3.1:

- `FVGridGeometry` and `EnableFVGridGeometryCache` have been replaced by
  `GridGeometry` and `EnableGridGeometryCache`.

Unfortunately, clang doesn't emit deprecation warnings if an old property name is
used. Consider employing gcc for detecting occurrences of the old name.

### Deprecated classes/files, to be removed after 3.1:

- `BaseFVGridGeometry` from `dumux/discretization`. Use `BaseGridGeometry` instead.

- `CakeGridCreator` and `SubgridGridCreator` from `dumux/io/grid`. Use the corresponding `...Manager` instead.

- `IntRange` from `dumux/common`. Use `Dune::IntegralRange` instead.

### Deprecated member functions, to be removed after 3.1:

- The convergence writer is no longer passed to `NewtonSolver`'s `solve()` method.
  For outputting convergence data, please use `newtonSolver.attachConvergenceWriter(convWriter)` in `main.cc` (directly after instantiating the writer).
  To stop the output, the writer can also be detached again using `newtonSolver.detachConvergenceWriter()`.

### Deleted classes/files, property names, constants/enums

- Deprecated classes and files for the 3.0 release listed below stay deprecated
  for at least one more release cycle.

Differences Between DuMu<sup>x</sup> 2.12 and DuMu<sup>x</sup> 3.0
=============================================

## Important Notes

- DuMu<sup>x</sup> 3.0 is a major version update. It is not backward compatible in all aspects to 2.12.
  The following minor version updates will be, as before for the DuMu<sup>x</sup> 2-series, always backward compatible
      to at least the last minor version update.

- The tutorial has been replaced by the new module `dumux-course` which is
  accessible at https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course.
  We recommend new users and also users experienced with DuMu<sup>x</sup> 2.X to clone
  the course and have a look at the exercises in there.

- DuMu<sup>x</sup> 3.0 is based on Dune 2.6 and is expected to run with the current Dune master.
  We will try to keep the compatibility with the Dune master
  as long as it is technically feasible and our resources allow it.

- DuMu<sup>x</sup> 3.0 requires at least GCC 4.9 or Clang 3.5 in their C++-14 mode.
  However, we suggest to use newer compiler versions, as we cannot test against all previous compiler versions.

- For employing corner-point grids by means of opm-grid, the OPM release 2018.10 has to be used.

## Improvements and Enhancements

### General

- __New style main files:__ 3.0 comes which a major overhaul of how the sequence of simulation steps is specified.
  Until 2.12 there was the start.hh header including the free function `start` that contained a predefined sequence of
  steps. Customization of the sequence was possible by many hooks implemented in classes used in this sequence. This
  made it hard to follow the course of the program and hard to implement customization at points where this was not previously envisioned.
  In contrast, in the new approach the sequence of simulation steps was linearized, meaning that each step is visible
  in the program's main file. In the new main files, one can clearly see the io, initialization, assembly, solving, update, time loop, etc.,
  of the program. Many parts of the simulation that were previously mandatory, i.e. simulations in 2.12 always contained a time loop,
  are now optional (stationary simulations just don't need a time loop). The new style main files make it easier two follow the course
  of the program, are easier to customize, and offer more flexibility concerning the customization of steps and components of a simulation.
  All tests and examples in the dumux repository have been adapted to the new style main files.
- __Property system:__ The property system is now usable without preprocessor macros. To this end it was completely reimplemented using C++14 techniques and
  variadic templates. The hierarchies can now be arbitrarily deep instead of being limited to 5 levels. The new implementation does not use
  C++ inheritance. Properties and TypeTag now have to be properly qualified with the namespaces `Properties::`, `TTag::`. Types that share the
  name with properties have to properly qualified with the `Dumux::` namespace. This update makes it hopefully more readable
  and removes the "magic" part from the property system.
- __Runtime parameters:__ Runtime parameters are no longer accessed with preprocessor macros. They have been replaced by C++ function templates
  `Dumux::getParam`, `Dumux::hasParam`, `Dumux::getParamFromGroup`. The `..FromGroup` version has been redesigned to allow the specification
  of parameters for different models in one input file. The concept of a parameter group string was extended to make it possible to
  use a single input file for complex multidomain simulation setups.
- __Restarting simulations:__ The old restart module was replaced by an implementation based on a VTK backend (other formats might be added later such as HDF5).
  Restarted simulations can read solutions from vtk files. For parallel runs, there is currently the restriction that the number of processors has to be the same
  as before the restart. Restarted simulations can read grids from vtk (currently only sequential, non-adaptive grids, support for parallel and adaptive
  will be added in a future version).
- __Assembly__: The assembler can now assemble implicit and explicit Euler time discretizations. An interface for implementing analytical Jacobians was added.
  The CCTpfa assembler has been significantly improved for complex models that spend a lot of time computing constitutive laws. Also the numerical
  differentiation scheme was improved by altering the order in which derivatives are computed.
- __TypeTag templates:__ Implementers of code in DuMu<sup>x</sup> 3.0 are advised to avoid TypeTag as a template argument for class templates.
  Many classes in the DuMu<sup>x</sup> core have been changed to have a small number of specific template arguments, including `GridGeometry`,
  `TimeLoop`, `Newton`, `LinearSolver`, `SpatialParams`. This makes it possible to share objects of these types between models using
  different TypeTags. This was not possible before as `Class<TypeTag1>` and `Class<TypeTag2>` are different types, even if they contain
  exactly the same implementation code.
  Furthermore, having TypeTag as a template argument leads to bad programming, and unnecessary dependencies that should be avoided in
  every object-oriented code.
- __The grid geometry concept:__ In version 3.0, discretization methods use grid geometries which are wrappers or adapters around a `Dune::GridView`,
  providing data structures and interfaces necessary to implement these discretization methods easily. In particular, the
  abstraction of sub-control-volumes (scv) and sub-control-volume-faces (scvf) are now separate classes and the grid geometry provides means to iterate
  over all scvs and scvfs of an element, using range-based-for loops.
- __The caching concept:__ Version 3.0 introduces a caching concept for the new grid geometry, the volume variables and the flux variables.
  There are classes with a `Grid` prefix, that store data for the current grid view, and classes with the `Element` prefix that store data locally
  for an element or element stencil. Depending on the caching concept used, data will be either stored completely in the `Grid` objects
  (e.g. `GridGeometry`) and the `Element` objects (e.g. `ElementGeometry`) are mere accessor objects, or, data will be partly only cached locally.
  Local caching uses less memory but might result in an increased runtime. Grid caching is memory intensive
  but can provide a significant run-time speedup. Choose whatever concept fits your available resources,
  the default being local caching.
- __Support__ for grid managers `dune-subgrid` (a meta grid selecting only certain elements from a host grid)
  and `dune-spgrid` (a structured parallel grid manager, supporting periodic boundaries).

### Models, Physics and Methods

- __MPFA schemes:__ The new design of the DuMu<sup>x</sup> core facilitates the incorporation of new finite-volume schemes. In particular, the new core comes with
  a framework for MPFA schemes, in which currently the only available scheme is the MPFA-O scheme. It can be used in conjunction with any DuMu<sup>x</sup> model and
  also works on surface grids. More schemes will be added in the future.
- __Box-dfm:__ The `2pdfm` model from version 2.12 has been generalized such that it can be used on any DuMu<sup>x</sup> model and in both two and three dimensions.
- __Tracer transport__: A new model for tracer transport with a given flow field has been added. The model can be also used to implement sequentially
  coupled simulations, or iterative solvers where flow and transport are decoupled / weakly coupled.
- __Mineralization__: An adapter model for mineralization has been added and can be used with all porousmediumflow models. A balance for the solid
  volume fraction of precipitating, adsorbed, or absorbed substances is added to the existing equations.
- __Solution-dependent spatial params:__ A redesign of the spatial params interface allows now to define spatial parameters such as permeability
  and porosity that depend on the solution. This makes it easier to implement mineralization models altering the solid structure of the porous medium.
- __Solid systems:__ DuMu<sup>x</sup> 3.0 introduces solid systems similar to fluid systems but for solid components. This allows a consistent
  implementation of mineralization models including reactions, dissolution, precipitation and other processes altering the solid
  phase of the porous medium.
- __Multidomain:__ DuMu<sup>x</sup> 3.0 introduces a new multidomain framework which does no longer depend on `dune-multidomain` and can be used for the coupling
  of an arbitrary number of subdomains. The sub-domains can be regions in which different sets of equations are solved and/or which have different
  dimensionalities. The implementation is such that any of the existing DuMu<sup>x</sup> models can be used in the subdomains, while the data and functionality
  required for the coupling of the sub-domains is implemented in a `CouplingManger` class. Three different coupling concepts are available, for which
  there are a number of available `CouplingManager` class implementations:
    - _Boundary:_ coupling across sub-domain boundaries
    - _Embedded:_ Coupling between a bulk domain and an embedded lower-dimensional sub-domain which has an independent grid
    - _Facet:_ Coupling between a bulk domain and a codimension-one sub-domain, which is conforming with the element facets of the bulk domain
- __Free-flow models:__ The previous Navier-Stokes model using the box method has been replaced by one that employs a staggered grid discretization.
  The new method does not  require any stabilization techniques while those were necessary for the box method in order to avoid spurious oscillations.
  The free-flow models in DuMu<sup>x</sup> 3.0 consider single phase flow with or without component/energy transport. So far, only regular Cartesian grids are supported
  but local grid refinement will be added in a future release.
  Several RANS models for turbulent flow have been added: k-omega, k-epsilon, low-Re-k-epsilon, one-eq, zero-eq. The RANS models might be subject to further (interface)
  changes.
- __Thermal and chemical non-equilibrium:__ The possibility to consider thermal and/or chemical non-equilibrium of several types has been enabled for all
  porous medium models.
- __Interface solver:__ For the two-phase flow model in conjunction with the box scheme, an interface solver can now be used to reconstruct the saturations
  in the sub-control volumes adjacent to vertices that lie on material discontinuities. This allows a sharper representation of the saturation front evolving
  in heterogeneous porous media.
- __Components:__ Components can now derive from different base classes, `Base`, `Liquid`, `Solid`, `Gas`, depending on which
  phase states are implemented. This can be used to determine at compile time if a component supports a certain phase state.
- __Tabulation of fluid parameter laws:__ The tabulation of fluid parameter laws has been improved to only tabulate those functions actually used during the
  simulation. To this end, the tabulation is done on the first call of a corresponding fluid parameter.
- __MPNC:__ The general m-phase n-component model has been adapted in structure to the other porous medium flow models.
- __Different wettability:__ The 2p models can now model materials with different wettability (hydrophobic, hydrophilic) in different parts of the domain.
- __Maxwell-Stefan-diffusion:__ Most models can now use Maxwell-Stefan diffusion for multi-component diffusion instead of Fickian diffusion.
  There are also a few tests demonstrating how to use it.

## Immediate interface changes not allowing/requiring a deprecation period

- Many classes have been completely redesigned. See the numerous example applications included in 3.0 showcasing all new classes.
- The `GridCreator` has been replaced by the `GridManager`, which no longer uses a singleton for the grid object.
  This makes it possible to create two grids of the exact same type. The `GridManager` also handles user data provided in grid files.

## Deprecated classes/files, to be removed after 3.0
- All classes of the sequential models are deprecated. The sequential models will be ported to the new structure
  of porous medium models (formerly called implicit models). This way sequential and implicit model implementations
  no longer differ and use the same numerical infrastructure.
- The `TimeManager` class is to be replaced by the class `TimeLoop`.
- The `VtkMultiWriter` class is to be replaced by the class `VtkOutputModule`.
- The file `start.hh` is replaced by new style main files.

## Deleted classes/files/functions... which have been deprecated in DuMu<sup>x</sup> 2.12
- Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.11 and DuMu<sup>x</sup> 2.12
=============================================

* IMPORTANT NOTES:
    - DuMu<sup>x</sup> 2.12 is expected to run based on Dune 2.4.1, 2.5, 2.6 and the Dune
      master. We will try to keep the compatibility with the Dune master
      as long as it is technically feasible and our resources allow it. If
      you want to use DuMu<sup>x</sup> multidomain models, you have to stick with the
      Dune 2.4 core and specific versions of other modules, see
      [test/multidomain/README](test/multidomain/README) for details.
      Also the geomechanics models require Dune 2.4 and PDELab 2.0.

    - DuMu<sup>x</sup> 2.12 requires at least GCC 4.9 or Clang 3.5 in their C++-14 mode.

    - For employing corner-point grids by means of opm-grid (former
      dune-cornerpoint), the OPM releases 2017.04 or 2017.10 have to be used.

* IMPROVEMENTS and ENHANCEMENTS:
    - Four new tutorial exercises have been added in the folder
      [tutorial](tutorial). They can be built by executing
      `make build_tutorials` in the build folder. Each exercise
      comes with detailed instructions:
        1. [Exercise 1](tutorial/ex1/README.md)
        2. [Exercise 2](tutorial/ex2/README.md)
        3. [Exercise 3](tutorial/ex3/README.md)
        4. [Exercise 4](tutorial/ex4/README.md)

    - Fixed bug in heatCapacity() of component air and replace
      the use of a constant value in  gasEnthalpy() by calling
      heatCapacity().

    - The GnuplotInterface now supports in-simulation generation of image
      files (*.png).

    - A paraview python script for exporting 2d pictures from *.vtu files
      has been added.

    - A class for estimating turbulence properties has been added with
      turbulenceproperties.hh.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
    - gnuplotinterface.hh: The add...ToPlot() functions have changed signature,
      the curve name/title is not mandatory anymore and can be specified together
      with the curve options.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.11:
    - Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.10 and DuMu<sup>x</sup> 2.11
=============================================

* IMPORTANT NOTES:
    - DuMu<sup>x</sup> 2.11 is expected to run based on Dune 2.4.1, 2.5 and the Dune
      master. We will try to keep the compatibility with the Dune master
      as long as it is technically feasible and our resources allow it. If
      you want to use DuMu<sup>x</sup> multidomain models, you have to stick with the
      Dune 2.4 core and specific versions of other modules, see
      `test/multidomain/README` for details.

    - DuMu<sup>x</sup> 2.11 requires at least GCC 4.9 or Clang 3.5 in their C++-14 mode.

    - For employing corner-point grids by means of opm-grid (former
      dune-cornerpoint), the OPM release 2016.04 has to be used.

* IMPROVEMENTS and ENHANCEMENTS:

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
    - A gridcreator for piece-of-cake-type grids has been added. It is capable
      of creating meshes with gradually in- and decreasing distances between nodes.
      It also allows the creation of a 360° cake where the last elements are
      connected to the first elements.
    - shouldWriteRestartFile() is now, as shouldWriteOutput() already was,
      called before the time level is advanced. So it might be necessary to use
      ...WillBeFinished instead of ...IsFinished for writing restart files at
      the correct time.
    - In the ZeroEq models, the properties BBoxMinIsWall and BBoxMaxIsWall have
      been replaced by the functions bBoxMaxIsWall() and bBoxMaxIsWall() in the
      problem file.
    - In the TwoPNC (and, consequently the TwoPNCMin) models, the old formulations
      pgSl, plSg as well as pnSw and pwSg have been replaced by the pnsw and pwsn,
      to satisfy the naming convention and be consistent with TwoPTwoC.
    - In the TwoPTwoC model, the indices are no longer dependent on the
      formulation. Further, the values of "nPhaseOnly" and "bothPhases"
      have been harmonized with those in TwoPNC

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.11: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. If possible, a run-time warning will appear in the summary lines
  after the corresponding run.

* Deprecated CLASSES/FILES, to be removed after 2.11:

* Deprecated MEMBER FUNCTIONS, to be removed after 2.11:

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.10:
    - Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.9 and DuMu<sup>x</sup> 2.10
===================================================

* IMPORTANT NOTES:
    - DuMu<sup>x</sup> 2.10 is expected to run based on Dune 2.4.1, 2.5 and the Dune
      master. We will try to keep the compatibility with the Dune master
      as long as it is technically feasible and our resources allow it. If
      you want to use DuMu<sup>x</sup> multidomain models, you have to stick with the
      Dune 2.4 core and specific versions of other modules, see
      `test/multidomain/README` for details.

    - DuMu<sup>x</sup> 2.10 requires at least GCC 4.9 or Clang 3.5 in their C++-14 mode.

    - For employing corner-point grids by means of opm-grid (former
      dune-cornerpoint), the OPM release 2016.04 has to be used.

* IMPROVEMENTS and ENHANCEMENTS:
    - Two new fully-implicit models have been added: `ThreePWaterOil` and
      `TwoPOneC`. The first one admits two components water and oil that may be
      present in two liquid and one gaseous phase, see `test_box3pwateroil` in
      `test/porousmediumflow/3pwateroil/implicit` for details. The second one is
      dedicated to one generic component that may be present in liquid and
      gaseous state, see `test_boxsteaminjection` in
      `test/porousmediumflow/2p1c/implicit` for an example.

    - Numbering of the VTK files starts with `0` now.

    - Using the geostatistical tool gstat for generating random fields has been
      facilitated. See `test_cc1pwithgstat` in `test/porousmediumflow/1p/implicit`.
      This tool can be installed using the `bin/installexternal.sh` script.
      If cmake does not find gstat, one has to specify the GSTAT_ROOT variable,
      see the standard optim.opts or debug.opts.

    - The multidomain models should now run with all compilers without
      segfaults, both with optimization and debug options.

    - Computation for two-dimensional fracture networks in three-dimensional
      space has been fixed and is tested now in
      `test/porousmediumflow/2p/implicit/test_fracture_box2p`.

    - The `bin` folder containing various helper scripts has been restructured.
      For example, the scripts `fuzzycompare.py` and `runtest.py` are now
      contained in the subfolder `testing`.

    - Two new scripts `extractlinedata.py` and `extractpointdataovertime.py`
      have been added to `bin/postprocessing`. They realize what their names
      suggest.

    - The computation of the tortuosity tau with the Millington-Quirk model has
      been optimized. Tests show a speedup of up to 5% for a 2p2c simulation.

    - The fully-implicit box scheme now fully supports prism-shaped elements.

    - Convenience functions for reading values from a file to a container and
      for writing values to a file from a container have been added, see
      `test_container_io` in `test/io/container`.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
    - The problem class interface needs to provide the method maxTimeStepSize().
      This is important if you implement problem classes not deriving from the base
      problem classes in DuMu<sup>x</sup> (`ImplicitProblem`, `OneModelProblem`,
      `ImpetProblem`, and `MultidomainProblem`).
    - All name-related methods that previously returned / received `const char*`
      do now use the type-safe alternative `std::string`. An example is
      `FluidSystem::componentName()`. If you need a
      `const char*` for special operation use the string member `c_str()`.

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.10: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. If possible, a run-time warning will appear in the summary lines
  after the corresponding run.
    - The pure run-time parameter `tau` has been renamed to the property and
      run-time parameter `TauTortuosity`. For the models 1p2c, 2p2c, 3p3c it got
      the default value 0.5.

* Deprecated CLASSES/FILES, to be removed after 2.10:
    - The class `DiffusivityConstantTau<TypeTag, Scalar>` was renamed to
      `DiffusivityConstant<TypeTag>` and is located in a according new header. The
      old class, its header and the property `tau` are deprecated.

    - `PointSourceHelper<TypeTag>` has been deprecated. Use
      `BoundingBoxTreePointSourceHelper<TypeTag>` instead.

    - The classes `Dumux::Liquid/GasPhase<Scalar, Component>` have been moved to
      the namespace `Dumux::FluidSystems`.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.10:
    - `checkEffectiveSaturation` of the classes `PlotMaterialLaw<TypeTag>` in
      `dumux/io`.

    - `concentrationGrad` of the class `TwoPNCFluxVariables<TypeTag>`. Use
      `moleFractionGrad` instead, with the obviously different physical meaning.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.9:
    - Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.8 and DuMu<sup>x</sup> 2.9
===================================================

* IMPORTANT NOTES:
    - DuMu<sup>x</sup> 2.9 is expected to run based on either Dune 2.4 or Dune 3.0. We will
      try to keep the compatibility with Dune 3.0 as long as it is technically
      feasible and our resources allow it. If you want to use DuMu<sup>x</sup> multidomain
      models, you have to stick with the Dune 2.4 core and specific versions of
      other modules, see `test/multidomain/README` for details.

    - The ALUGrid stand-alone library cannot be used any longer. Use the module
      [dune-alugrid](https://gitlab.dune-project.org/extensions/dune-alugrid)
      instead. Both the `releases/2.4` branch and the `master` should work.

    - The above is not true if you like to run a sequential MPFA model on an
      `ALUGrid`. Then, you currently have to use the `master-old` branch of
      dune-alugrid. We will try to fix this as soon as possible. Alternatively,
      `UGGrid` or `YaspGrid` can be chosen as grid managers.

    - Instead of using AlbertaGrid for the tests where dim < dimWorld, we now
      employ
      [dune-foamgrid](https://gitlab.dune-project.org/extensions/dune-foamgrid).
      Dune-foamgrid provides 1d and 2d simplex grids embedded in an arbitrary
      dimension world space. It features element parametrizations, runtime growth,
      runtime-movable vertices. You might still use AlbertaGrid, but it is not
      supported by our GridCreator.

    - If you like/have to use corner-point grids by means of the module
      dune-cornerpoint, you have to use (and partially patch) the 2015.10 release
      of [OPM](http://opm-project.org/?page_id=36). See `patches/README` for
      details.

* IMPROVEMENTS and ENHANCEMENTS:
    - The folder structure has been changed according to
      [FS#250](http://www.dumux.org/flyspray/index.php?do=details&task_id=250).
      This has been a rather massive change affecting more than 1000 files. Close
      to 400 files have been moved and/or renamed.
      We made everything backwards-compatible, the worst thing that should happen
      after switching to DuMu<sup>x</sup> 2.9, will be some warnings when including headers
      from old destinations/names. You can fix the include statements and get rid
      of the warnings by applying the bash script `bin/fix_includes.sh` to your
      source files, for example by executing
      ```
      bash ../dumux/bin/fix_includes.sh file1 [file2 ...]
      ```
      or
      ```
      find . -name '*.[ch][ch]' -exec bash ../dumux/bin/fix_includes.sh {} \;
      ```
      inside the folder that contains your files.
      A patch is available to remove deprecated header files:
      ```
      patch -p1 < patches/dumux-2.9-no-deprecated-headers.patch
      ```
      The benefits are hopefully:
        + A clearer structure in terms of the problems that you want to apply DuMu<sup>x</sup>
          for. Three main application areas on the top level: `porousmediumflow`,
          `freeflow` and `geomechanics`. The different numerical treatments "fully
          implicit" or "sequential" appear as discretization detail after the
          choice of the physical model. That's of course currently rather wishful
          thinking, but nevertheless where we are headed. The folder `implicit` on
          the top level now only contains physics-agnostic classes that can be used
          by any class of an application area.  Please note the change from
          "decoupled" to "sequential" according to the related task
          [FS#252](http://www.dumux.org/flyspray/index.php?do=details&task_id=252).

        + Nicer include statements due to relaxation of the naming conventions for
          the header files. Compare the old
          ```
          #include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnilocaloperator.hh>
          ```
          with the new
          ```
          #include <dumux/multidomain/2cnistokes2p2cni/localoperator.hh>
          ```
      The structure change is reflected in the `test` folder:
        + The tests from`test/implicit/particular_model` have been moved to
          `test/porousmediumflow/particular_model/implicit`. For example,
          `test/implicit/2p` has been moved to `test/porousmediumflow/2p/implicit`.

        + Analogously, the tests from `test/decoupled/particular_model` have been
          moved to `test/porousmediumflow/particular_model/sequential`.

        + The subfolders `decoupled` and `implicit` of `test` have been removed.

        + If you have cloned the DuMu<sup>x</sup> Git repository and have local changes in the
          folders `test/implicit` or `test/decoupled`, you can expect merge
          conflicts for your next `git pull`. You can either deal with these
          conflicts directly or create a patch, remove the local changes, pull, and
          apply the patch afterwards with some care to respect the changed
          structure.

    - A two-phase multiple-interacting-continua (MINC) model has been added to
      the DuMu<sup>x</sup> model portfolio. See `test/porousmediumflow/2pminc/implicit` for
      details.

    - The multidomain models have been restructured. Duplicated code has been
      reduced; isothermal and non-isothermal models are treated in a more
      consistent manner.

    - It is now possible to specify point sources for implicit models. A point
      source is a source term specified at any point location in e.g. kg/s. DuMu<sup>x</sup>
      will compute the correct control volume the source belongs to for you. Point
      sources can be e.g. solution and/or time-dependent. See tests
      (1p/implicit/pointsources, 2p/implicit/pointsources) for examples.

    - All tests use our standard `GridCreator` now. If it is possible to specify
      the grid entirely in the input-file, the corresponding DGF files have been
      deleted. In particular, a YaspGrid tensor grid can now also be specified via
      the input file only.

    - Several sections on our fluid/material framework have been moved from the
      handbook to the Doxygen documentation.

    - The three-phase constitutive relations from
      `material/fluidmatrixinteractions` have been reworked to be consistent with
      their two-phase analogues. In particular, an `EffToAbsLaw` and
      regularization classes have been implemented.

    - In case of a simulation stop due to too many timestep subdivisions, restart
      files of both the current and the old solution are automatically generated.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
    - All flux variables are now default-constructible. While the non-trivial
      constructors are deprecated, model implementers might be required to make
      their flux variables default-constructible too. In particular, this affects
      you if you develop your own flux variables that
        + inherit from flux variables from dumux-stable, such as the
          `ImplicitDaryFluxVariables`,
        + and/or are used in a local residual from dumux-stable.
      See the
      [mailing list](https://listserv.uni-stuttgart.de/pipermail/dumux/2016q1/001551.html)
      for details.

    - For the multidomain models, the notation of the boundary condition types
      has changed. This is especially important for all momentum boundary
      conditions. In general:
        + `couplingInflow`  -> `couplingNeumann`
        + `couplingOutflow` -> `couplingDirichlet`
    - But for the momentum balances:
        + `couplingInflow`  -> `couplingDirichlet`
        + `couplingOutflow` -> `couplingNeumann`

    - Due to the change in the three-phase fluid-matrix-interactions, you might
      have to adjust your spatial parameters. You should get a compiler warning
      message that gives you more details.

    - The TypeTags `ImplicitModel` and `ExplicitModel` have been deleted. They
      haven't been used apart from one internal inheritance. See FS#304 for
      details.

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.9: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.
    - The word `Decoupled` in the TypeTags has been replaced by `Sequential`:
        + `DecoupledModel` -> `SequentialModel`
        + `DecoupledOneP` -> `SequentialOneP`
        + `DecoupledTwoP` -> `SequentialTwoP`
        + `DecoupledTwoPTwoC` -> `SequentialTwoPTwoC`
        + `DecoupledTwoPTwoCAdaptive` -> `SequentialTwoPTwoCAdaptive`


* Deprecated CLASSES/FILES, to be removed after 2.9:
    - Self-written parallel linear solvers and corresponding infrastructure,
      according to
      [FS#293](http://www.dumux.org/flyspray/index.php?do=details&task_id=293).
      For parallel runs, use the `AMGBackend` instead. For sequential runs,
      direct replacements are:
        + `BoxBiCGStabILU0Solver` -> `ILU0BiCGSTABBackend`
        + `BoxBiCGStabSORSolver` -> `SORBiCGSTABBackend`
        + `BoxBiCGStabSSORSolver` -> `SSORBiCGSTABBackend`
        + `BoxBiCGStabJacSolver` -> `JacBiCGSTABBackend`
        + `BoxBiCGStabGSSolver` -> `GSBiCGSTABBackend`
        + `BoxCGILU0Solver` -> `ILUnCGBackend`
        + `BoxCGSORSolver` -> `SORCGBackend`
        + `BoxCGSSORSolver` -> `SSORCGBackend`
        + `BoxCGJacSolver` -> `JacCGBackend`
        + `BoxCGGSSolver` -> `GSCGBackend`
        + `IMPETBiCGStabILU0Solver` -> `ILU0BiCGSTABBackend`

    - `CubeGridCreator`, functionality available in default `GridCreator`

    - `SimplexGridCreator`, functionality available in default `GridCreator`

    - `DgfGridCreator`, functionality available in default `GridCreator` (since 2.8)

    - `Decoupled...Indices` -> `Sequential...Indices` (BEWARE: maybe no compiler
    warnings)

* Deprecated MEMBER FUNCTIONS, to be removed after 2.9:

* Deprecated protected MEMBER VARIABLES, to be removed after 2.9: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable
  is used.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.8:
    - Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.7 and DuMu<sup>x</sup> 2.8
===================================================

* IMPORTANT NOTES:
  - DuMu<sup>x</sup> 2.8 is expected to run based on either Dune 2.3 or Dune 2.4. However,
    no attempt has been made to fix the warnings arising from the deprecation of
    EntityPointer in Dune 2.4. This will be made during the next release cycle.
    Moreover, using the multidomain models based on Dune 2.4 is currently only
    possible by patching dune-multidomaingrid. See test/multidomain/README for
    details.

* DELETED BUILD SYSTEM: The Autotools based build system was removed, use the
  CMake based build system as it is default since Dune 2.4.

* IMPROVEMENTS and ENHANCEMENTS:
  - New fully-implicit porous-media models for two fluid phases that may consist
    of an arbitrary number of components have been added. The basic one is
    associated with the property TwoPNC, see test/implicit/2pnc. A more advanced
    one that incorporates solid-fluid phase changes is indicated by TwoPNCMin,
    see test/implicit/2pncmin.

  - The implicit cell-centered models now can use adaptive grid refinement. To
    make a test problem adaptive, just set the property AdaptiveGrid to true and
    choose corresponding indicators via AdaptionInitializationIndicator and
    AdaptionIndicator, see test/implicit/2p/lensproblem.hh for an example. So
    far, indicators are only provided for the TwoPModel. Indicators for other
    models will be provided in the future, as well as parallelization and box
    discretization.

  - With the CpGridCreator, a grid creator has been introduced that reads from a
    Petrel output / Eclipse input file and generates a CpGrid that is offered by
    the OPM module dune-cornerpoint. The fully-implicit cell-centered models are
    now able to deal with cornerpoint grids. See
    test/implicit/2p/test_cc2pcornerpoint for a test of the functionality. A
    realistic corner-point grid will be provided in dumux-lecture soon. The OPM
    modules need to be patched to be compatible with Dune's CMake based build
    system, see patches/README for details.

  - Zero equation turbulence models (zeroeq) have been added as new models
    to the freeflow folder. Tests for coupling a turbulent free flow using
    zeroeq turbulence models with flow in a porous medium have been added.

  - A new class GridCreator is now the new standard grid creator replacing
    DgfGridCreator. It comprises all functionality from the DgfGridCreator (see
    also immediate interface changes), plus the ability to read gmsh, or to
    build a structured grid (only with Dune 2.4) by merely changing the input
    file.

  - Multidomain problems can now be run by using the general start routine,
    just as most other problems. For this, the constructor of the multidomain
    problems has been changed and the InterfaceMeshCreator has been replaced by
    the InterfaceGridCreator, see below.

  - The Richards model has now an additional flag useHead, which can be used to
    switch between a pressure-saturation and a pressureHead-watercontent
    formulation. The primary variables are either pressure in [Pa] or pressure
    head in [cm], respectively. Default is useHead = false. See
    test/implicit/richards for details.

  - A bug in the diffusion term in the freeflow models has been fixed.

  - A lot of work has been devoted to improving the testing environment, adding
    new tests, restructuring the handbook and improving the documentation.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:

  - The (new standard) GridCreator's public method gridPtr has been removed. For
    using PARAMETERS from DGF files use the GridCreators parameters method via
    GridCreator::parameters(...). See test/implicit/co2 for an example.

  - The use and support for SGrid is dropped. SGrid is deprecated in Dune 2.4.
    Use YaspGrid instead.

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.8: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.

* Deprecated CLASSES/FILES, to be removed after 2.8:
  - The InterfaceMeshCreator has been moved to InterfaceGridCreator and adapted
    to the structure of other grid creators, it can simply be used by specifying
    the GridCreator TypeTag.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.8:
  - The constructor of the multidomain problems has changed. They will now be
    called directly from the start.hh, identical to the other problems. The new
    parameters are the TimeManager and the HostGrid.

  - The method simulate(Scalar dtInitial, Scalar tEnd) from MultiDomainProblem,
    is unused and will be dropped.

  - The GnuplotInterface functions are now called without giving a window
    number. If plots should be plotted in different windows, different
    GnuplotInterface objects are now required. This affects also all other plots
    in the "io" folder.

  - The write() function in plotoverline2d.hh now has an append function, to be
    able to decide whether the previously written file should be kept or not.

* Deprecated protected MEMBER VARIABLES, to be removed after 2.8: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable
  is used.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.7:
  Everything listed as deprecated below has been removed.

Differences Between DuMu<sup>x</sup> 2.6 and DuMu<sup>x</sup> 2.7
===================================================

* IMPORTANT NOTES:
  - DuMu<sup>x</sup> 2.7 should work with DUNE 2.3 as well as 2.4. However, at the time of
    writing, DUNE-multidomain(grid) doesn't work with DUNE 2.4. Therefore, if a
    DuMu<sup>x</sup> multidomain model should be used, DUNE 2.3 is required. See
    test/multidomain/README for details.

  - The 2.3 branch of dune-alugrid has no CMake support, use dune-alugrid master
    respectively 2.4. Or you can fall back to Autotools or use legacy ALUGrid
    1.52.

* IMPROVEMENTS and ENHANCEMENTS:
  - Since 2.6, all isothermal implicit porous-media models (except 2pdfm) can be
    easily enhanced with the energy equation. For 2.7, this has been also
    carried out for the models that were only isothermal before, namely, 1p, 3p
    and richards. Tests have been written and are provided in test/implicit. In
    order to keep the number of subfolders bearable, isothermal as well as
    thermal tests are gathered in the model folders "1p", "1p2c", ..., "3p3c",
    "mpnc", "richards" (without the "ni") and the corresponding "ni"-folders
    have been deleted.

  - All implicit porous-media models (except 2pdfm) are now able to run on grids
    with dim < dimWorld. In implicit/1p, four new tests are added that run the
    1p test problem on 1d-3d and 2d-3d Alberta grids with box and cell-centered,
    respectively. Compilation has been tested also for all other models, but no
    runtime testing has been performed.

  - The terminology for the Newton method has been improved according to FS#238.
    In particular, what has been referred to as "relative error" is now termed
    "maximum relative shift", while "absolute error" has been renamed to
    "residual reduction". This is particularly important, if corresponding
    parameters or properties are set, see below.

  - The geomechanics ElTwoPModel runs in parallel now. This is made possible by
    a dedicated solver, the El2PAMGBackend which has to be set for the property
    LinearSolver in the problem file. See test/geomechanics/el2p for details.

  - Before, velocity output for the implicit porous-media models only worked for
    cube grids. This has been generalized to simplices (box and cc) and prisms/
    pyramids (box only).

  - Revised and fixed restart capability for the multidomain models.

  - A gnuplot interface has been added (works only with CMake). With
    this interface it is possible to plot material laws (like in the 2p2cni test),
    or to generate live-updating output (like in test_2cnistokes2p2cni).
    The gnuplot interface reads analytical functions, data file or data arrays.

  - The fuzzycompare script for automatic testing has been improved. Instead of
    printing only the first deviation from the reference solution, it now
    prints the maximum deviation in each field/variable.

* DEPRECATED BUILD SYSTEM: DuMu<sup>x</sup> 2.7 will be the last release which can be built
    with the Autotools based build system. It is deprecated and will be removed
    for DuMu<sup>x</sup> 2.8. We encourage the change towards CMake, especially with the
    upcoming DUNE 2.4.
    The warning can be suppressed with --disable-dumux-deprecated-autotools

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - Before, the "heatCapacity" function in the spatial parameters and volume
    variables of the implicit nonisothermal models was a misnomer, since it
    returned an effective quantity, namely,
    heatCapacity*density*(1 - porosity) in [J/(K m^3)].
    Except for mpnc, which resulted in an additional inconsistency.
    Corresponding to the decision documented in FS#216, the function has been
    renamed to "solidHeatCapacity" and returns always the "true" (non-effective)
    heat capacity in [J/(kg K)]. This requires an additional function
    "solidDensity" which returns the mass density of the porous matrix.
    Moreover, the functions "thermalConductivitySolid/Fluid" are renamed to
    "solid/fluidThermalConductivity". The decision to prepend with "solid/fluid"
    rather than to append is motivated by consistency with components and fluid
    systems, where "gas" and "liquid" are always prepended to the corresponding
    function names.
    Therefore, it might be necessary to adapt your thermal solid parameters in
    the spatialparams file such that they offer functions "solidHeatCapacity",
    "solidDensity" and "solidThermalConductivity". See
    test/implicit/2p2c/injectionspatialparams.hh for an example.

  - Due to the change in the Newton terminology (see above), there exist two
    backward-compatibility breakages:
    . If a model re-implements the function "relativeErrorDof", it has to be
      renamed to "relativeShiftAtDof". See dumux/implicit/implicitmodel.hh for
      an example.

    . If a NewtonController re-implements the function "newtonUpdateRelError",
      it has to be renamed to "newtonUpdateShift". See
      dumux/nonlinear/newtoncontroller.hh for an example.

  - The properties "AMGPDELabBackend" and "AMGLocalFemMap" have been unified to
    "AmgTraits".

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.7: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.
  - Corresponding to the improved Newton terminology, the following properties
    (prepended with "Newton") and parameters (in the group "Newton") are
    renamed:
    AbsTolerance -> ResidualReduction
    EnableAbsoluteCriterion -> EnableResidualCriterion
    RelTolerance -> MaxRelativeShift
    EnableRelativeCriterion -> EnableShiftCriterion
    SatisfyAbsAndRel -> SatisfyResidualAndShiftCriterion

* Deprecated CLASSES/FILES, to be removed after 2.7:
  - SeqAMGBackend and ScaledSeqAMGBackend, replaced by AMGBackend.

  - P0LocalFiniteElementMap.

  - CellData2P2Cmultiphysics, replaced by CellData2P2CMultiPhysics.

  - BoxLocalOperator from dumux/multidomain/common/pdelablocaloperator.hh.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.7:
  - The functions "heatCapacity", "densitySolid" (mpnc only) and
    "thermalConductivitySolid/Fluid" in the VolumeVariables of the nonisothermal
    implicit porous-media models: use "solidHeatCapacity", "solidDensity" and
    "solid/fluidThermalConductivity" instead. See also the immediate interface
    changes above.

  - In dumux/implicit/common/implicitmodel.hh and
    dumux/geomechanics/el2p/elp2basemodel.hh:
    "relativeErrorDof" -> "relativeShiftAtDof"

  - In dumux/nonlinear/newtoncontroller.hh:
    "setRelTolerance" -> "setMaxRelativeShift"
    "setAbsTolerance" -> "setResidualReduction"
    "newtonUpdateRelError" -> "newtonUpdateShift"

  - The 1p2c volume variables no longer use the method tortuosity() from
    spatial params class, the value is now calculated within the effective
    diffusivity model. Thus the method is deprecated in the spatial params
    classes FVSpatialParamsOneP and ImplicitSpatialParamsOneP.

* Deprecated protected MEMBER VARIABLES, to be removed after 2.7: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable
  is used.
  - In dumux/nonlinear/newtoncontroller.hh:
    "error_" -> "shift_"
    "lastError_" -> "lastShift_"
    "tolerance_" -> "shiftTolerance_"
    "absoluteError_" -> "reduction_"
    "lastAbsoluteError_" -> "lastReduction_"
    "initialAbsoluteError_" -> "initialResidual_"
    "absoluteTolerance_" -> "reductionTolerance_"
    "enableRelativeCriterion_" -> "enableShiftCriterion_"
    "enableAbsoluteCriterion_" -> "enableResidualCriterion_"
    "satisfyAbsAndRel_" -> "satisfyResidualAndShiftCriterion_"

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.6:
  Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.5 and DuMu<sup>x</sup> 2.6
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - For the non-isothermal porous media models, the energy equation is
    implemented in a more generic way for all models in
    dumux/implicit/nonisothermal. The existing TypeTag names like TwoPNI stay
    the same. If a new non-isothermal model should be used, it is important
    to NOT include anything from the old model-specific implementation like
    from dumux/implicit/2pni, but to include from the model folder without the
    "ni". See test/implicit/2pni for details. In principle, any isothermal
    porous media model can be enhanced with the energy equation. Ideally, only
    the corresponding property files have to be augmented. See
    dumux/implicit/2p/2ppropert*.hh for details. The 1p2c model already has
    been enhanced, the remaining models will follow in 2.7.

  - The AMG backend is based directly on dune-istl now. No PDELab is required
    anymore. The tests so far exhibit an improved robustness. Thanks to Markus
    Blatt for the work.

  - The multidomain models can now be used with the 2.3 release versions of the
    DUNE core modules and dune-multidomaingrid, and the 2.0 release versions
    of dune-pdelab and dune-multidomain. See test/multidomain/README for
    details.

  - In the fully implicit mpnc model, a further specialization allows now to
    describe two-phase flow with two energy equations.

  - The free flow models now include the component enthalpy fluxes transported
    by diffusion processes (h^k D grad x), which was not considered before.

  - UMFPack is a new direct linear solver and can be use as a drop-in
    replacement for SuperLU. Some users claim a speed-up up to a factor of
    seven. We know cases where it was 10% slower, so please measure for your
    problems.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - The 3p3cni model now also uses an effective thermal conductivity model
    (ETCM). The ETCM is easily exchangeable. The default one is
    ThermalConductivitySomerton, which is implemented in
    dumux/material/fluidmatrixinteractions/3p. The ETCM requires that 3p3cni
    spatial parameters provide a function thermalConductivitySolid instead of
    matrixHeatFlux. See test/implicit/3p3cni/columnxylolspatialparams.hh for
    details. Moreover, the employed fluid system has to actually implement the
    function thermalConductivity. See
    dumux/material/fluidsystems/h2oairxylenefluidsystem.hh for details.

  - The non-isothermal flux variables call the effective thermal conductivity
    models (ETCM) in a different way. If you used a self-written ETCM and want
    to use a new non-isothermal model, the ETCM has to be adapted. See
    material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh for
    details.

  - Fully implicit mpnc model: in order to account for the possibility of using
    two energy equations, the boolean property EnableKineticEnergy has been
    changed to the integer property NumEnergyEquations.

* Deprecated way of setting command line parameters, to be removed after 2.6:
  - To set parameters from the command line, the notation --parameterFile=NAME
    is deprecated. Use from now on -ParameterFile NAME.

* Deprecated CLASSES/FILES, to be removed after 2.6:
  - The old non-isothermal porous media models are deprecated. Technically,
    including a ..niproperties.hh file triggers a deprecation warning.

  - FVPressure2P2CAdaptive, use dimension-specific implementations
    FV2dPressure2P2CAdaptive and FV3dPressure2P2CAdaptive instead.

  - FVTransport2P2CAdaptive, use dimension-specific implementations
    FV2dTransport2P2CAdaptive and FV3dTransport2P2CAdaptive instead.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.6:
  - In the Stokes flux variables, the method eddyViscosity() is deprecated, use
    dynamicEddyViscosity() instead.

  - In the Stokes non-isothermal flux variables, the method eddyConductivity()
    is deprecated, use thermalEddyConductivity() instead.

  - Already in 2.5, the following member functions of MultiDomainModel/Problem
    have been deprecated: subProblemX, subModelX, subIDX, gridViewX with X=1,2.
    They are replaced by sdProblemX, sdModelX, sdIDX, sdGridViewX.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.5:
  Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.4 and DuMu<sup>x</sup> 2.5
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - The three-dimensional implementation of the MPFA L-method is made
    available for the decoupled compositional 2p2c models. It also allows
    for simulation on an adaptive grid.
  - Coupling of 2c2p with stokesnc and 2p2cNI with Stokesncni was
    added. The stokes2c and stokes2cni are now DEPRECATED and will be kicked
    out by the next release. Instead generalized stokesnc and stokesncni
    models are introduced. Unlike 2c models the transport equations in
    the nc models are capapable of using both mass and mole fractions.
    NOTE: For coupling test examples be aware of the harsh version
    restrictions mentioned in dumux/test/modelcoupling/README.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - Dropped support for PDELab 1.0.1, PDELab 1.1 is required.

* Deprecated CLASSES/FILES, to be removed after 2.5:
  - Stokes2cModel was replaced by StokesNCModel, similar for more
    Stokes2c* classes.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.4:
  Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.3 and DuMu<sup>x</sup> 2.4
===================================================

* IMPORTANT NOTES:
  - If the current trunk version of DUNE (2.3) is used, the co2 and co2ni tests
    require that the DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS flag is set. This
    flag is needed in order to be able to use the function boundaryId(). An
    error is thrown during compilation if the flag is not set. For reference see
    the commented part in optim.opts and debug.opts.

  - All two-component models (1p2c, 2p2c, 2p2cni, co2, co2ni) can now be used
    with either mole or mass fractions. The property useMoles has to be set in
    the problem file and the boundary conditions have to be chosen accordingly.
    . 1p2c, 2p2c, 2p2cni use MOLE fractions by default.
    . co2, co2ni use MASS fractions by default.
    . For completeness: 3p3c, 3p3cni and mpnc use MOLE fractions only.

* IMPROVEMENTS and ENHANCEMENTS:
  - Three geomechanics models have been added, thanks to the previous PhD work
    of Melanie Darcis. The models are located in dumux/geomechanics, a test for
    each model is provided in test/geomechanics:
    . elastic: linear elasticity equations that account only for the solid
               mechanics.
    . el1p2c:  poroelasticity equations for one solid and one fluid phase, where
               the fluid phase is admitted to consist of two components.
    . el2p:    poroelasticity equations for one solid phase and two immiscible
               fluid phases.

  - A three-dimensional implementation of the MPFA L-method is added and made
    available for the decoupled 2p models. It also allows for simulation on
    an adaptive grid.

  - The MPNC model now allows to drop the assumptions of local thermal and/or
    chemical equilibrium. Dropping the local chemical equilibrium assumption
    leads to mole fractions in different phases that are not determined via
    equilibrium relations. If the local thermal equilibrium assumption is
    dropped, phases (fluids and solid) are locally allowed to have different
    temperatures. See test/implicit/mpnc/test_boxmpnckinetic for an example.

  - A fully-implicit three-phase immiscible model has been added. The headers
    are contained in dumux/implicit/3p, tests are provided in test/implicit/3p.

  - The handling of Dirichlet boundary conditions for the fully-implicit cell-
    centered models has been improved. Now, the conditions are evaluated at the
    centers of the corresponding boundary faces. In general, they are
    incorporated into the local residual in a weak sense. Only for mixed
    Dirichlet/Neumann conditions and for the MPNC model, they still are
    incorporated in a strong sense.

  - The sequential models can use a sub-timestepping in the transport scheme, if
    the property "SubCFLFactor" is chosen smaller than "CFLFactor", which in
    that case specifies the CFL factor used in the pressure equation.

  - All fully-implicit porous-media models now provide the possibility to write
    the velocities to the VTK output. This can be achieved by setting the
    parameter "AddVelocity" of the group "Vtk" or the corresponding property
    "VtkAddVelocity." See test/implicit/1p2c for an example.

  - The CMake build system support uses the experimental mechanisms from DUNE.
    Recent versions of Dune 2.3-svn or newer are required to use CMake.

  - Our naming rules have been refined and enforced for member functions and
    variables. Many inconsistencies could be removed, especially the special
    status of the capitalized "S" indicating saturation. See the deprecation
    listings below or FS#180 for details.

  - Misleading names in the fully-implicit models that still contained "b/Box"
    have been replaced by more generic ones. See the deprecation listings below
    or FS#194 for details.

  - The ...FVElementGeometry classes have been cleaned up a bit. See the
    deprecation listings below or FS#192 for details.

  - Added compiler support for Clang 3.2, Clang 3.3, and GCC 4.8.

* Deprecated CLASSES/FILES, to be removed after 2.4:
  - OnePBoxModel, OnePTwoCBoxModel -> OnePModel, OnePTwoCModel

  - Headers containing capitalized letters are replaced by their non-capitalized
    analogs.
    In dumux/decoupled: 1p/cellData1p.hh, 1p/fluxData1p.hh,
    2p/cellData2padaptive.hh, 2p/fluxData2p.hh, 2p/cellData2p.hh,
    2p2c/cellData2p2c.hh, 2p2c/cellData2p2cadaptive.hh, 2p2c/fluxData2p2c.hh,
    2p2c/cellData2p2cmultiphysics.hh.
    In dumux/material/fluidmatrixinteractions/3p: parkerVanGen3p.hh,
    parkerVanGen3pparams.hh.

* Deprecated CONSTANTS/ENUMS, to be removed after 2.4: BEWARE: Older compilers
  will not print any warning if a deprecated constant/enum is used.
  - saturation indices: S0Idx, SnIdx, SnOrXIdx, SOrX1Idx, SOrX2Idx, SwIdx,
    SwOrXIdx
    -> s0Idx, snIdx, snOrXIdx, sOrX1Idx, sOrX2Idx, swIdx, swOrXIdx

  - two-phase formulations: pnSn, pnSw, pwSn, pwSw -> pnsn, pnsw, pwsn, pwsw

  - DecoupledTwoPCommonIndices: pressureNW, saturationNW, velocityNW
    -> pressureNw, saturationNw, velocityNw

  - DecoupledTwoPIndices: pressEqIdx -> pressureEqIdx

  - MPNCIndices: NumPrimaryEnergyVars, NumPrimaryVars
    -> numPrimaryEnergyVars, numPrimaryVars

* Deprecated public MEMBER VARIABLES, to be removed after 2.4: BEWARE: Older
  compilers will not print any warning if a deprecated public member variable
  is used.
  - ...FVElementGeometry: numEdges, numFaces, numFap, numVertices
    -> numScvf, -, numFap of each subcontrolvolume face, numScv

  - BoxFVElementGeometry: edgeCoord, faceCoord

* Deprecated MACROS, to be removed after 2.4: BEWARE: The compiler will not
  print any warning if a deprecated macro is used.
  - DUMUX_ALWAYS_INLINE whether the according attribute is supported

* Deprecated MEMBER FUNCTIONS, to be removed after 2.4:
  - all problems: bboxMin/Max() -> bBoxMin/Max()

  - ImplicitProblem: boxSDNeumann(), boxSDSource()
    -> solDependentNeumann(), solDependentSource()

  - ImplicitPorousMediaProblem: boxGravity(), boxTemperature()
    -> gravityAtPos(), temperatureAtPos() (different signatures!)

  - fluid-matrix-interactions: dkrn_dSw(), dkrw_dSw(), dpc_dSw(), pC(),
    dSw_dpC(), Sgr(), Snr(), SnToSne(), Sw(), Swr(), SwToSwe()
    -> dkrn_dsw(), dkrw_dsw(), dpc_dsw(), pc(), dsw_dpc(), sgr(), snr(),
    snToSne(), sw(), swr(), swToSwe()

  - LinearMaterial(Params): entryPC(), maxPC(), setEntryPC(), setMaxPC()
    -> entryPc(), maxPc(), setEntryPc(), setMaxPc()

  - RegularizedVanGenuchten(Params): pCHighSw(), pCLowSw()
    -> pcHighSw(), pcLowSw()

  - VanGenuchtenParams, ParkerVanGen3PParams: setVgM(), setVgN(), vgM(), vgN()
    -> setVgm(), setVgn(), vgm(), vgn()

  - ParkerVanGen3P(Params): betaGN(), betaGW(), betaNW(), pCAlpha(), pCGN(),
    pCGW(), pCNW(), setBeta..., setkrRegardsSnr(), Swrx()
    -> betaGn(), betaGw(), betaNw(), pcAlpha(), pcgn(), pcgw(), pcnw(),
    setBeta..., setKrRegardsSnr(), swrx()

  - MPLinearMaterialParams: Sreg() -> sReg()

  - EvalCflFlux...: getCFLFluxFunction() -> getCflFluxFunction

  - FVMPFAOInteractionVolume: getNTKNu_by_dF(), getNTKNu(), getNTKrKNu_by_dF(),
    getNTKrKNu()
    -> getNtkNu_df(), getNtkNu(), getNtkrkNu_df(), getNtkrkNu()

  - TwoPDFMVolumeVariables: dSM_dSF() -> dsm_dsf()

  - Stokes...Variables: viscosity() -> dynamicViscosity()

  - IAPWS water: ddgamma_ddpi, ddgamma_ddtau, ddgamma_dtaudpi, dgamma_dpi,
    dgamma_dtau, dp_dpi, dpi_dp, dtau_dt
    -> ddGamma_ddPi, ddGamma_ddTau, ddGamma_dTaudPi, dGamma_dPi, dGamma_dTau,
    dp_dPi, dPi_dp, dTau_dt

* DELETED classes/files, property names, constants/enums,
  member functions/variables, which have been deprecated in DuMu<sup>x</sup> 2.3:
  Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.2 and DuMu<sup>x</sup> 2.3
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - A fully implicit two-phase discrete-fracture-matrix model has been added,
    see test/implicit/2pdfm.

  - Almost all porous media fully implicit models now can either use a
    vertex-centered (box) or a cell-centered spatial discretization. The choice
    of the spatial discretization method is controlled by deriving the problem
    type tag either from BoxModel or CCModel. This allows for a uniform problem
    description, as long as the boundaryTypesAtPos and dirichletAtPos methods
    can be used. By evaluating the compile-time property ImplicitIsBox, it is
    easily possible to separately handle the different discretizations inside
    am common method. See the tests in test/implicit for examples.
    Correspondingly, the directory structure has been adapted.
    Old:             New:
    dumux/           dumux/
      boxmodels/       implicit/
        common/          common/
        1p/              box/
        1p2c/            cellcentered/
        2p/              1p/
        ...              ...
    test/            test/
      boxmodels/       implicit/
        1p/              1p/
          test_1p          test_box1p
        ...                test_cc1p
                         ...

  - A backend for the ISTL AMG solver has been included, based on the
    corresponding DUNE-PDELab backends. It can be used for the fully
    implicit and the decoupled models, see test_*1pwithamg in
    test/implicit/1p and test_impeswithamg in test/decoupled/2p.
    DUNE-PDELab and possibly DUNE-ISTL have to be patched, see the file
    README in the patches directory.

  - The decoupled models have been parallelized, see test_impeswithamg in
    test/decoupled/2p. They work in parallel only if the AMGBackend is used
    as linear solver. No dynamic loadbalancing can be done yet.

  - The MPNC model can use either the most wetting or the most nonwetting phase
    pressure as primary variable. This is controlled via the property
    "PressureFormulation."

  - The table of available parameters has been improved, see
    http://www.dumux.org/doxygen-stable/html-2.2/a00838.php

  - Improved handling of the conductive heat fluxes in the non-isothermal implicit
    two-phase models, see the problem files in test/implicit/2p(2c)ni.

  - Introduced new selection of start/stop messages.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - The property Salinity used in the BrineCO2FluidSystem
    has been renamed to ProblemSalinity.

  - The matrixHeatFlux(...) and boundaryMatrixHeatFlux(...) methods in the spatial
    parameters of nonisothermal implicit twophase models have been removed.
    Instead, the computation of the effective thermal conductivity has been sourced
    out to the fluidmatrixinteractions in a separate file
    dumux/material/fluidmatrixinteractions/thermalconductivitysomerton.hh, which
    can be exchanged. The spatial parameters file needs a method
    thermalConductivitySolid(...), where the thermal conductivity of the solid
    material only is specified. The rest is computed in the respective
    flux variables.

* Deprecated CLASSES/FILES, to be removed after 2.3:
  - The following headers in dumux/boxmodels/ have been deprecated and forward
    to the corresponding headers in dumux/implicit/box:
    boxassembler.hh, boxelementvolumevariables.hh, boxlocalresidual.hh,
    boxpropertydefaults.hh, boxelementboundarytypes.hh, boxfvelementgeometry.hh,
    boxproperties.hh, intersectiontovertexbc.hh

  - All headers in the following subdirectories of dumux/boxmodels have been
    deprecated and forward to the headers in the corresponding subdirectories
    of dumux/implicit:
    1p, 1p2c, 2p, 2p2c, 2p2cni, 2pdfm, 2pni,
    3p3c, 3p3cni, co2, co2ni, mpnc, richards

  - Some box-specific classes "Box..." in dumux/boxmodels/common could be
    completely replaced by unified "Implicit..." classes in
    dumux/implicit/common:
    ...DarcyFluxVariables, ...darcyfluxvariables.hh
    ...ForchheimerFluxVariables, ...forchheimerfluxvariables.hh
    ...LocalJacobian, ...localjacobian.hh
    ...Model, ...model.hh
    ...PorousMediaProblem, ...porousmediaproblem.hh
    ...Problem, ...problem.hh
    ...VolumeVariables, ...volumevariables.hh

  - The box-specific spatial parameter classes BoxSpatialParams... in
    dumux/material/boxspatialparams....hh have been deprecated in favor of
    ImplicitSpatialParams... in dumux/material/implicitspatialparams....hh.

  - The GridCreatorheaders from dumux/common have been moved to dumux/io:
    cubegridcreator.hh, dgfgridcreator.hh, simplexgridcreator.hh

* Deprecated PROPERTY NAMES, to be removed after 2.3: BEWARE: The compiler will
  not print any warning if a deprecated property name is used.
  - CompositionFromFugacitiesSolver has been renamed to Constraintsolver.

* Deprecated public MEMBER VARIABLES, to be removed after 2.3: BEWARE: The
  compiler will not print any warning if a deprecated public member variable
  is used.
  - numFAP and numSCV in Box(CC)FVElementGeometry have been renamed to
    numFap and numScv, respectively.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.3:
  - boundaryMatrixHeatFlux, markVertexRed and relativeErrorVertex
    from ImplicitSpatialParams, ImplicitAssembler and ImplicitModel,
    respectively. In favor of using
    thermalConductivitySolid (see above), markDofRed and relativeErrorDof,
    respectively.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMu<sup>x</sup> 2.2:
  Everything listed as deprecated below has been removed.


Differences Between DuMu<sup>x</sup> 2.1 and DuMu<sup>x</sup> 2.2
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - Two new fully implicit models dedicated to simulate compositional
    (non-isothermal) CO2-brine systems have been added, together with
    corresponding components and a fluid system. See test/boxmodels/co2(ni)
    for details. These tests also illustrate the usage of element and vertex
    parameters as well as boundary ids provided by DGF files for setting
    permeability and porosity as well as boundary conditions.

  - Decoupled Models: An h-adaptive model using an MPFA L-method was added
    that simulates 2p and 2p2c flow on unstructured grids with hanging nodes
    in two dimensions. See test/decoupled/2p(2c) for details.

  - All fully implicit porous media models are now capable of employing
    the Forchheimer equation as an alternative to the commonly used
    Darcy law. See test_forchheimer*p in test/boxmodels/mpnc for details.

  - The Stokes models are now able to simulate the full Navier-Stokes
    equations for momentum transport. See test/freeflow/navierstokes
    for details.

  - The fully implicit models have been (partially) generalized to allow
    for a cell-centered discretization in addition to the default
    vertex-centered (box) one. Cell-centered fully implicit 2p and 2p2c
    models are already available in the developers part of DuMu<sup>x</sup>. Further
    generalizations and the inclusion in the stable part are planned for
    DuMu<sup>x</sup> 2.3.

  - Several model-specific features and classes have been unified, like
    the calculation of the Darcy velocity for the fully implicit flux
    variables, or the temperature, gravity, and spatial parameter
    functionalities of the fully implicit problems. Moreover, many
    names have been made more consistent. This includes the naming
    and grouping of several parameters and corresponding properties,
    the indexing of phases and components, and the preference of the
    partial name "params" over "parameters." For details, see also the
    deprecations listed below.

  - Added compiler support for GCC 4.7 and Clang 3.1.

* IMMEDIATE INTERFACE CHANGES not allowing a deprecation period:
  - From Dune 2.2 on, FieldVector::size is a method rather than an enum value.
    It is mandatory to add the flag --enable-fieldvector-size-is-method to the
    CONFIGURE_FLAGS. An example is given in the opts file dumux/debug.opts.
  - Implicit models: TwoPIndices, TwoPNIIndices, and RichardsIndices
    additionally get TypeTag as template parameter. If the Indices are not
    obtained via the property, this has to be adapted.

  - Implicit models: All model-specific computeFlux functions in
    ...localresidual.hh have to get an additional bool parameter onBoundary,
    which is by default set to false. If outflow conditions should
    be properly implemented, also the constructor of the flux variables in
    ...fluxvariables.hh has to get the additional argument and the
    class has to be adapted to deal with boundary faces. See FS#117 and #99
    for details.

* Deprecated CLASSES/FILES, to be removed after 2.2:
  - Model specific base box problems: The common functionality has been
    collected in PorousMediaBoxProblem in
    dumux/boxmodels/common/porousmediaboxproblem.hh. The problem can be derived
    from PorousMediaBoxProblem, instead of the model specific base problem:
    OnePBoxProblem, dumux/boxmodels/1p/1pproblem.hh,
    OnePTwoCBoxProblem, dumux/boxmodels/1p2c/1p2cproblem.hh,
    TwoPProblem, dumux/boxmodels/2p/2pproblem.hh,
    TwoPNIProblem, dumux/boxmodels/2pni/2pniproblem.hh,
    TwoPTwoCProblem, dumux/boxmodels/2p2c/2p2cproblem.hh,
    TwoPTwoCNIProblem, dumux/boxmodels/2p2cni/2p2cniproblem.hh,
    ThreePThreeCProblem, dumux/boxmodels/3p3c/3p3cproblem.hh,
    ThreePThreeCNIProblem, dumux/boxmodels/3p3cni/3p3cniproblem.hh,
    MPNCProblem, dumux/boxmodels/mpnc/mpncproblem.hh.

  - All "...SpatialParameters" base classes have been replaced by
    "...SpatialParams" classes:
    BoxSpatialParameters, dumux/material/spatialparameters/boxspatialparameters.hh,
    BoxSpatialParametersOneP, dumux/material/spatialparameters/boxspatialparameters1p.hh,
    FVSpatialParameters, dumux/material/spatialparameters/fvspatialparameters.hh,
    FVSpatialParametersOneP, dumux/material/spatialparameters/fvspatialparameters1p.hh.

  - Due to the unification of flux variables for the fully implicit models,
    some model-specific flux variables have become obsolete:
    OnePFluxVariables, dumux/boxmodels/1p/1pfluxvariables.hh,
    TwoPFluxVariables, dumux/boxmodels/2p/2pfluxvariables.hh,
    RichardsFluxVariables, dumux/boxmodels/richards/richardsfluxvariables.hh.

  - Two components have new names and locations in dumux/material/components:
    SimpleDNAPL, simplednapl.hh -> DNAPL, napl.hh
    Oil, oil.hh -> LNAPL, lnapl.hh

  - Some MPFA-O method files/classes have been moved to a new subdirectory
    "omethod" in dumux/decoupled/2p/diffusion/fvmpfa:
    fvmpfaopressure2p.hh, fvmpfaovelocity2p.hh, fvmpfaopressureproperties2p.hh

  - DUMUX_UNUSED is deprecated and will be removed after 2.2. It should be
    replaced by the upstream version DUNE_UNUSED.

  - DUMUX_DEPRECATED_MSG is deprecated and will be removed after 2.2. It should
    be replaced by the upstream version DUNE_DEPRECATED_MSG.

* Deprecated PROPERTY NAMES, to be removed after 2.2: BEWARE: The compiler will
  not print any warning if a deprecated property name is used.
  - The "SpatialParameters" property has been renamed to "SpatialParams".

  - The model specific "...Indices" property has been renamed to "Indices".

* Deprecated CONSTANTS/ENUMS, to be removed after 2.2: BEWARE: The compiler will
  not print any warning if a deprecated constant/enum is used.
  - In the 2p2c/ni and 3p3c/ni models, all indices related to phase and
    components can be pre/suffixed with "w", "n" and,
    for three phases, with "g".
    boxmodels/2p2c/...: "l", "g" pre/suffixes have been replaced by "w", "n".
    boxmodels/3p3c/...: "c", "a" pre/suffixes have been replaced by "n", "g".

* Deprecated MEMBER FUNCTIONS, to be removed after 2.2:
  - Spatial parameters: The spatialParameters member functions of the base
    problems have been replaced by spatialParams:
    dumux/boxmodels/common/porousmediaboxproblem.hh,
    dumux/decoupled/1p/diffusion/diffusionproblem...hh,
    dumux/decoupled/2p/impes/impesproblem2p.hh,
    dumux/decoupled/2p/transport/transportproblem2p.hh.

  - Flux variables: Renaming of members
    "...AtIP" -> "...",
    "concentration..." -> "massFraction...",
    "molarConc..." -> "moleFraction..."
    The "massFraction..." members have been deprecated, instead
    "moleFraction..." should be used.
    Affected files:
    dumux/boxmodels/1p2c/1p2cfluxvariables.hh,
    dumux/boxmodels/2p2c/2p2cfluxvariables.hh,
    dumux/boxmodels/mpnc/.../...fluxvariables.hh,
    dumux/freeflow/stokes.../stokes...fluxvariables.hh.

  - Box models: The primaryVarWeight() functions are no longer used for the
    evaluation of the relative error.

  - Element and FVElementGeometry: The elem_() and fvElemGeom_() member function
    of BoxLocalResidual have been replaced by element_() and fvGeometry_().

  - Primary variables: All "...primaryVar/s" member functions have been replaced
    by "...priVar/s":
    dumux/boxmodels/common/boxlocalresidual.hh,
    dumux/boxmodels/common/boxvolumevariables.hh.

  - Start functionality in dumux/common/start.hh: printUsageDGF and
    printUsageGrid are no longer needed.

* DELETED member functions, which have been deprecated in DuMu<sup>x</sup> 2.1:
  - dumux/material/spatialparameters/boxspatialparameters1p.hh:
    extrusionFactorScv and extrusionFactorScvf, now part of the volume variables

  - dumux/material/idealgas.hh:
    concentration, replaced by molarDensity

  - dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh:
    pC(const Params &params, Scalar Sw, const Scalar temperature)

  - dumux/common/start.hh:
    startFromDGF, startWithGrid, startWithParameters, all replaced by start

  - dumux/common/spline.hh:
    set(const ScalarArray&, const ScalarArray&, Scalar, Scalar), replaced by setXYArrays,
    set(const PointArray&, Scalar, Scalar), replaced by setArrayOfPoints

  - dumux/common/variablelengthspline_.hh, dumux/common/fixedlengthspline_.hh:
    various set routines, replaced by more descriptive names

  - dumux/io/vtkmultiwriter.hh:
    VtkMultiWriter(const std::string&, std::string), replaced by VtkMultiWriter(const GridView&, ...),
    beginTimeStep, replaced by beginWrite,
    createField, replaced by allocateManagedBuffer,
    addVertexData, replaced by attachVertexData,
    addCellData, replaced by attachCellData,
    endTimeStep, replaced by endWrite

  - dumux/decoupled/2p2c/2p2cproblem.hh:
    IMPETProblem2P2C(const GridView&, bool) replaced by IMPETProblem2P2C(TimeManager&, ...),
    IMPETProblem2P2C(..., SpatialParameters&, ...) replaced by IMPETProblem2P2C(TimeManager&, ...),
    initSat(const GlobalPosition&, const Element&) replaced by initSat(const Element&)
    initConcentration(const GlobalPosition&, const Element&) replaced by initConcentration(const Element&)

  - DUMUX_DEPRECATED has been removed.


Notable Differences Between DuMu<sup>x</sup> 2.0 and DuMu<sup>x</sup> 2.1
===================================================

- The thermodynamics framework has been overhauled:
  - The programming interfaces for fluid systems, fluid states and
    components has been formalized and cleaned up.
  - Fluid systems now have the option to cache computationally
    expensive parameters if they are needed for several relations.
  - Fluid systems are not charged with the computation of the
    chemical equilibrium anymore.
  - Fluid states are now centralized infrastructure instead of being
    model-specific.
  - Constraint solvers (which simplify solving thermodynamic
    constraints) have been introduced.
- Outflow boundary conditions have been implemented for the
  fully-implicit models 1p2c, 2p2c(ni) and stokes(2cni).
- The problem and spatial parameter base classes also provide optional
  model-independent interfaces. These methods only get the position in
  global coordinates as argument and are named *AtPos()
  (e.g. boundaryTypesAtPos()). This allows an easy transfer of problem
  definitions between implicit and sequential models.
- The following fully-implicit models have been added:
  - 3p3c, 3p3cni: Isothermal and non-isothermal three-phase,
    three-component models for flow and transport in porous media
    based on primary variable switching.
  - MpNc: A model for arbitrary number of phases M > 0, and components
    (N >= M - 1 >= 1) for flow and transport in porous media. This
    model also comes with an energy and a molecular diffusion module.
  - stokes, stokes2c, stokes2cni: Models for the plain Stokes
    equation as well as isothermal and non-isothermal Stokes models
    for two-component fluids.
- The sequentially-coupled models have been overhauled:
  - A common structure for cell centered standard finite volume
    implementations has been introduced.
  - The data structures where overhauled to avoid large clumps of data
    in large-scale simulations: Each cell stores data in its own
    storage object.
  - The too large assemble() methods have been split into submethods
    getStorage(), getFlux() etc. By this, inheritance of classes has
    been improved and code duplication was reduced.
  - The conceptual separation of the "VariableClass" (central
    infrastructure), data storage, transport model and pressure model
    has been improved.
  - More of infrastructure is now shared with the implicit models
    (e.g. the BoundaryTypes). This results in significant performance
    improvements and makes debugging easier.
- The 2padaptive sequentially coupled model has been added. This model
  implements a grid-adaptive finite volume scheme for immiscible
  two-phase flow in porous media on non-conforming quadrilateral
  grids.
- The dependencies for the external dune-pdelab and boost packages
  have been removed.
- The build system has received major improvements:
  - There is now much better test coverage of build-time dependencies
    on packages for the default autotools-based build system.
  - Experimental support for building DuMu<sup>x</sup> using CMake has been much
    improved. In the long run, CMake is projected to become the
    default build system.
- All headers can now be included without any preconditions.
- DuMu<sup>x</sup> now compiles without warnings if the -pedantic flag used for GCC.
- Specifying run-time parameters is now possible. The mechanism allows
  to use parameter files or to specify parameters directly on the
  command line and fallback parameter input files have been added for
  each test application.  As a consequence, applications can be run
  now without specifying any command line arguments.
- The DuMu<sup>x</sup> property system has been fine-tuned:
  - Encapsulating property names with the PTAG() is no longer required
    for the GET_PROP* macros (but is still allowed).
  - Setting property defaults has been deprecated.
  - All properties defined for a type tag can now be printed. Also,
    their value and the location in the source where they where
    specified is included in the output.
- Using quadruple precision math has been made possible for GCC 4.6 or newer:
  - To use it, the configure option '--enable-quad' needs to be added
    and the type of scalar values needs to be changed to 'quad'. This
    can be done in the problem file using

        SET_TYPE_PROP(YourProblemTypeTag, Scalar, quad);

    It should be noted, that performance is very poor when using
    quadruple precision arithmetic. This feature is primarily meant as
    a debugging tool to quickly check whether there are machine
    precision related convergence problems.
