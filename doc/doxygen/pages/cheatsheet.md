# Cheatsheet

## C++ Code
* Iterating over grid elements:
```cpp
for (const auto& element : elements(gridGeometry.gridView()))
{
    auto elementPos = element.geometry().center();
    ...
}
```

* Iterating over vertices:
```cpp
for (const auto& vertex : vertices(gridGeometry.gridView()))
{
    auto vertexPos = vertex.geometry().corner(0);
    ...
}
```

* Iterating over local vertices of an element:
```cpp
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
constexpr int dim = GridView::dimension;
for (const auto& element : elements(gridGeometry.gridView()))
{
    for (const auto& vertex : Dune::subEntities(element, Dune::Codim<dim>{}))
    {
        const auto vertexIdx = gridGeometry.vertexMapper().index(vertex);
        ...
    }

}
```

* Frequent calls on an element:
```cpp
element.level();                    // refinement level
auto geometry = element.geometry();
geometry.center();                  // center position
geometry.volume();                  // volume

// access global element index
auto elementIdx = gridGeometry.elementMapper().index(element);
// access element given global element index
auto element = gridGeometry.element(elementIdx);
```

* Frequent calls on a grid geometry:
```cpp
gridGeometry.numScv();                      // total number of sub-control volumes
gridGeometry.numScvf();                     // total number of sub-control volume faces
```
A grid geometry can also be accessed through a problem object:
```cpp
problem.gridGeometry.numScv();
```

* Iterating over all sub-control volumes:
```cpp
for (const auto& element : elements(gridGeometry.gridView()))
{
    const auto fvGeometry = localView(gridGeometry).bind(element);
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
scv.indexInElement(); // element-local index of the scv
```

* Creating and accessing volume variables:
```cpp
auto fvGeometry = localView(gridGeometry);
auto elemVolVars = localView(gridVariables->curGridVolVars());
for (const auto& element : elements(gridGeometry.gridView()))
{
    fvGeometry.bindElement(element);
    elemVolVars.bindElement(element, fvGeometry, solutionVector);
    for (const auto& scv : scvs(fvGeometry))
    {
        const auto& volVars = elemVolVars[scv];
        Scalar pressure = volVars.pressure(phaseIdx);
        Scalar density = volVars.density(phaseIdx);
    }
}

```
`gridVariables` is usually defined in the main file

* Creating and accessing flux volume variables, example for CCTPFA and porous medium flow model:
```cpp
using FluxVariables =  GetPropType<TypeTag, Properties::FluxVariables>;
auto upwindTerm = [](const auto& volVars) { return volVars.mobility(Indices::conti0EqIdx); };
auto fvGeometry = localView(gridGeometry);
auto elemVolVars = localView(gridVariables.curGridVolVars());
auto elemFluxVars = localView(gridVariables.gridFluxVarsCache());
for (const auto& element : elements(gridGeometry.gridView()))
{
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, solutionVector);
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

* Adding a data field to the output writer. In the following, the output writer `vtkWriter` is assumed to be of type `Dumux::VtkOutputModule` or `Dumux::IO::OutputModule`. The data field container needs to outlive the vtkwriter object:
```cpp
std::vector<double> outputField(gridGeometry.numDofs(), 0.0);
vtkWriter.addField(outputField, "f");
// The writer stores a reference: outputField cannot go out of scope
// outputField can be updated later, without updating the writer
```

* Adding a specific volume variable, e.g. the porosity, to the output writer. `volVars` is the argument of a lambda function and does not need to be defined before this call:
```cpp
vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.porosity(); }, "porosity");
```

## Working with (DUNE) modules
* Configure and compile all modules (including DuMux) using specified CMake options:
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
```
* Build only the module `<module_dir>`:
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --only=<module_dir> all
```
* Out-of-source build in directory `<absolute_build_directory_path>`:
```bash
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts --build-dir=<absolute_build_directory_path> all
```

* Removing cache before reconfiguring:
```bash
./dune-common/bin/dunecontrol bexec rm -rf CMakeFiles CMakeCache.txt
```

* Running with MPI (distributed-memory parallelism) and configuring shared-memory parallelism:
```bash
DUMUX_NUM_THREADS=<num_threads_per_core> mpirun -np <n_cores> ./<executable_name>
```

* Only configuring shared-memory parallelism:
```bash
DUMUX_NUM_THREADS=<num_threads> ./<executable_name>
```

* Passing run-time parameters per command line:
```bash
./dummy_test -GroupName.ParameterName parameterValue
./dummy_test <parameter_file> -GroupName.ParameterName parameterValue
```

* Build documentation locally. In `build-cmake/`:
```bash
make doc
```
Entry page of built documentation can be found in `dumux/build-cmake/doc/doxygen/html/index.html` and opened via a browser.

Build and run the test suite of DuMux. In `build-cmake/`:
```bash
make -j8 build_test & ctest -j8
```

* Creating a new DUNE module:
```bash
./dune-common/bin/duneproject
```

* Extract test(s) as extra module:
```bash
./dumux/bin/extractmodule/extract_as_new_module.py <module_dir> <sub_folder_1> <[sub_folder_2 ...]>
```

* Create install script for a module. Is automatically executed within `extract_as_new_module.py`:
```bash
./dumux/bin/extractmodule/make_installscript.py -p <module_dir>
```

* Create docker image of a module:
```bash
./dumux/bin/extractmodule/create_dockerimage.py -m <module_dir> -i <path_to_installscript>
```

## CMake (Build system)
* Setting debug option within the `CMakeLists.txt` file for a test, options for `[BUILD_TYPE]`: `Debug`, `RelWithDebInfo`:
```cmake
set(CMAKE_BUILD_TYPE [BUILD_TYPE])
```

## Useful Resources
- **DuMuX website:**
  [https://dumux.org/](https://dumux.org/)
- **DuMuX documentation:**
  [https://dumux.org/docs/doxygen/master/](https://dumux.org/docs/doxygen/master/)
- **DuMux style guide:**
  [https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/doc/styleguide.md](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/doc/styleguide.md)
- **Contribution guidelines:**
  [https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md)
- **DuMuX Wiki:**
  [https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/home](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/home)
- **DuMuX course:**
  [https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course)
- **Installation instructions:**
  [https://dumux.org/docs/doxygen/master/installation.html](https://dumux.org/docs/doxygen/master/installation.html)
- **List of all runtime parameters:**
  [https://dumux.org/docs/doxygen/master/parameterlist_8txt.html](https://dumux.org/docs/doxygen/master/parameterlist_8txt.html)
- **DUNE cheatsheet**:
  [https://www.dune-project.org/pdf/dune-cheat-sheet.pdf](https://www.dune-project.org/pdf/dune-cheat-sheet.pdf)

## General Remarks
- Clarify whether an object is a pointer or a reference before calling functions on it.
  - Use `->` for pointers and `.` for references when calling functions on objects.
  - Example: `gridGeometry->gridView()` vs. `gridGeometry.gridView()`
