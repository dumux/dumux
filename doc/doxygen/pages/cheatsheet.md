# Cheatsheet

## CMake Commands
* Setting debug option within a test, options for `[BUILD_TYPE]`: `'Debug'`, `'RelWithDebInfo'`:

```cmake
set(CMAKE_BUILD_TYPE [BUILD_TYPE])
```

## Command Line Commands
* Configuring modules:
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
```
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --only=<moduleName> all
```
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --build-dir=<buildDirectoryPath> all
```

* Removing cache before reconfiguring:
```bash
./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt
```

* Removing build folders before reconfiguring:
```bash
rm -rf d*/build-cmake
```

* Running executable with MPI (distributed-memory parallelism):
```bash
mpirun -np [n_cores] [executable_name]
```

* Running executable with shared-memory parallelism:
```bash
DUMUX_NUM_THREADS=[num_threads] ./test_executable
```

* Creating a new Dune module:
```bash
./dune-common/bin/duneproject
```

* Passing run-time parameters per command line:
```bash
./dummy_test -GroupName.ParameterName parameterValue
```

* Build documentation locally. In `build-cmake/`:
```bash
make doc
```
Entry page of built documentation can be found in `dumux/build-cmake/doc/doxygen/html/index.html` and opened via a browser.

## C++ Code
* Iterating over grid elements:
```cpp
for (const auto& element : elements(problem.gridGeometry().gridView()))
    // do something
```

* Frequent calls on an element:
```cpp
element.geometry().center();  // center position
element.geometry().volume();  // volume
element.level();              // refinement level

// access global element index
unsigned int elementIdx = gridGeometry.elementMapper().index(element);
```

* Frequent calls on a grid geometry:
```cpp
gridGeometry.numScv();   // total number of sub-control volumes
gridGeometry.numScvf();  // total number of sub-control volume faces
```

* Iterating over all sub-control volumes:
```cpp
for (const auto& element : Dune::elements(gridGeometry.gridView()))
{
    auto fvGeometry = localView(gridGeometry);
    fvGeometry.bind(element);

    for (const auto& scv : scvs(fvGeometry))
    {
        // do something
    }
}
```

* Frequent calls on a sub-control volume:
```cpp
scv.center();        // center position
scv.volume();        // volume
scv.dofIndex();      // index of DOF the scv is embedded in
scv.localDofIndex(); // element-local index of DOF the scv is embedded in
scv.dofPosition();   // position of DOF the scv is embedded in
```

* Creating and accessing volume variables:
```cpp
using VolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::VolumeVariables;

VolumeVariables volVars;
const auto elemSol = elementSolution(element, solutionVector, gridGeometry);
volVars.update(elemSol, problem, element, scv);
Scalar pressure_phase = volVars.pressure(phaseIdx);
Scalar density_phase = volVars.density(phaseIdx);
```

* Creating and accessing flux volume variables, CCTPFA example:
```cpp
using FluxVariables =  GetPropType<TypeTag, Properties::FluxVariables>;
auto upwindTerm = [](const auto& volVars) { return volVars.mobility(Indices::conti0EqIdx); };

auto fvGeometry = localView(gridGeometry);
auto elemVolVars = localView(gridVariables.curGridVolVars());
auto elemFluxVars = localView(gridVariables.gridFluxVarsCache());
for (const auto& element : elements(gridGeometry.gridView()))
{
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, sol);
    elemFluxVars.bind(element, fvGeometry, elemVolVars);

    for (const auto& scvf : scvfs(fvGeometry))
    {
        //treat inner scvfs
        if (!scvf.boundary())
        {
            FluxVariables fluxVars;
            fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
            auto innerFlux = fluxVars.advectiveFlux(Indices::conti0EqIdx, upwindTerm);
        }
        //treat boundary scvfs
        else
        {
            const auto bcTypes = boundaryTypes(element, scvf);
            if (bcTypes.hasOnlyDirichlet()) // Dirichlet
            {
                FluxVariables fluxVars;
                fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                auto boundaryFluxD = fluxVars.advectiveFlux(Indices::conti0EqIdx, upwindTerm);
            }

            else // Neumann
            {
                auto boundaryFluxNM = problem.neumann(element, fvGeometry, elemVolVars, 0.0, scvf)[Indices::pressureIdx]
                                          * scvf.area() * elemVolVars[0].extrusionFactor()
                                          / elemVolVars[0].density(); // volume flux from mass flux
            }
        }
    }
}
```

* Adding data fields to the output writer:
```cpp
vtkWriter.addField(<dataField>, "<fieldName>");
```

## Useful Resources
- **DuMuX website:**
  [https://dumux.org/](https://dumux.org/)
- **DuMuX documentation:**
  [https://dumux.org/docs/doxygen/master/](https://dumux.org/docs/doxygen/master/)
- **DuMuX course:**
  [https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course)
- **Installation instructions:**
  [https://dumux.org/docs/doxygen/master/installation.html](https://dumux.org/docs/doxygen/master/installation.html)
- **List of all runtime parameters:**
  [https://dumux.org/docs/doxygen/master/parameterlist_8txt.html](https://dumux.org/docs/doxygen/master/parameterlist_8txt.html)

## General Remarks
- Ensure that an object is a pointer or a reference before calling functions on it.
  - Use `->` for pointers and `.` for references when calling functions on objects.
  - Example: `gridGeometry->gridView()` vs. `gridGeometry.gridView()`

- Solution vector structure for Darcy model:
```cpp
{{pw1, sn1},
 {pw2, sn2},
 ...
 {pwN, snN}}
```