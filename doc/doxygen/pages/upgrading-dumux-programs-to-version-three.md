# Upgrading dumux programs to version three

The major version update from DuMuX 2.12 to DuMuX 3.0 included some backward-incompatible changes. For the 3-series backward-compatibility will again be assured from one minor version update to the next.

In the following guide we try to describe in detail which interfaces changed and how they changed from 2.12 to 3.0, concerning the `Problem` class, the `SpatialParams` class, the `Parameters` and `Properties`, i.e. the main user interface. Due to the large number of changes, we cannot discuss all changes in other classes in Dumux in this detail.

If there is anything missing, please report this on the mailing list (dumux@listserv.uni-stuttgart.de), or propose changes to this guide.

## Property system
The new property system can be used without macros. Most of the old macros are still available (except diagnostic macros). There is one difference in the usage: `NEW_PROP_TAG(...)` is now a _definition_ instead of a _declaration_ and can thus only appear once. You can remedy this problem by collecting all `NEW_PROP_TAG(...)` calls of your model/problem in a common header file with a header guard. You can also replace the first occurrence of macro by the actual property definition

```cpp
template<class TypeTag, class MyTypeTag>
struct NameOfMyProperty {
    //...default implementation... or
    //using type = Undefined;
};
```
and use forward declarations elsewhere. How to use the property system without macros is described in the handbook. Furthermore, all example applications in Dumux now use the new property system without macros. Just have a look at one of the `problem.hh` files.


## Parameters
The parameter tree is now implemented as a singleton.
Runtime parameter are now obtained with the free functions
``` cpp
template<class ParamType>
ParamType getParam<ParamType>(const std::string&);
```
and
``` cpp
template<class ParamType>
ParamType getParamFromGroup<ParamType>(const std::string&, const std::string&);
```
The functions optionally take a third argument providing a default if the parameter was not found in the parameter tree, i.e.
``` cpp
template<class ParamType>
ParamType getParam<ParamType>(const std::string&, ParamType defaultValue);
```

Furthermore, the tree can be queried for existence of a key with

``` cpp
bool hasParam(const std::string&);
```

The parameter tree is initialize by calling the `Parameters::init` function.
The `init` function exists with different signatures (see `dumux/common/parameters.hh`).
Most commonly it is called at program start like

``` cpp
// initialize the parameters
Parameters::init(argc, argv);
```

This default call looks first for a parameter file with the name of executable `<executablename>.input`, then for `params.input`, and finally it prints a message that no parameter file could be found and continues the program without a parameter file. Then, command line arguments are read, they overwrite the parameters in the parameter file, e.g.

``` bash
./executablename -TimeLoop.TEnd 1e6
```
sets or overwrites the parameter `"TimeLoop.TEnd"`.
Some parameters have global defaults which can be found in `dumux/common/parameters.hh`.
If a runtime parameter was never specified, a `ParameterException` is thrown.

## Main file
In DuMux 2.12 the main files typically only contained a call to `Dumux::start()`, a convenience function to carry out the standard program flow related to instationary, fully-implicit and non-linear problems (provided in the header `start.hh`). This is abandoned in DuMuX 3.0 in favor of writing out the major steps of the program flow in the main file, which provides more flexibility and better readability. Consequently, you are no longer restricted to use the framework for instationary, fully-implicit and non-linear problems, but it is now possible to solve e.g. stationary and linear problems without having to use neither a time loop nor a non-linear solver (see e.g. [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/test/porousmediumflow/1p/implicit/incompressible/main.cc)).

Particularly, this means that the main simulation flow is no longer hierarchically structured, and relying on hooks in different classes that are called in a predefined order. Instead, the full main simulation control is written in a linear fashion in the main file. This also means that many functions in the `Problem` class related to simulation flow (what step comes when) have been removed (see below).


## Problem class
The base classes for problems, originally `ImplicitPorousMediaProblem` and `ImplicitProblem`, are now independent of the time-stepping scheme used and are provided in the classes `PorousMediumFlowProblem` and `FVProblem` (finite volume problem). Being usable in both stationary and instationary simulations, the constructors no longer expect a `TimeManager` (which was replaced by the `TimeLoop` class in DuMuX 3.0), but solely require a `std::shared_ptr` to an `FVGridGeometry`. The latter is a wrapper around a `Dune::GridView` for the construction of the geometries of the chosen finite volume scheme, and consequently the following functions which originally were part of the __problem__ interface have been moved to this class:

`grid()`, `gridView()`, `bBoxMin()`, `bBoxMax()`, `vertexMapper()`, `elementMapper()`, `dofMapper()`

The `FVGridGeometry` and it's element-wise local counter part `FVElementGeometry` provide the abtractions of sub control volumes and sub control volume faces. For instance it it now possible given an instance of `FVElementGeometry` called `fvGeometry` to iterate over all sub control volumes associated with this element

```cpp
for (const auto& scv : scvs(fvGeometry))
{
    std::cout << scv.center() << std::endl;
}
```


With the new style of main files, allowing users to take action at any point along the program flow, the interfaces `preTimeStep()` and `postTimeStep()`, originally designed for intervention in certain points during the time integration, have been removed. The same holds for the interfaces related to grid adaption, i.e. `preAdapt()`, `postAdapt()` and `gridAdapt()`. Furthermore, __time control__, __solving of the linear/nonlinear system__, __restart mechanisms__ and __input/output__ now occur on the uppermost level within the `main` function, so that the following interfaces have been removed as well:

`timeIntegration()`, `newtonMethod()`, `newtonController()`, `nextTimeStepSize()`, `shouldWriteRestartFile()`, `shouldWriteOutput()`, `maxTimeStepSize()`, `advanceTimeLevel()`, `episodeEnd()`, `currentVTKFileNumber()`, `timeManager()`, `serialize()`, `restart()`, `deserialize()`, `addOutputVtkFields()`, `writeOutput()`, `resultWriter()`.

You can check out how the new __restart facility__ works by taking a look at e.g. [this](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/test/porousmediumflow/2p2c/implicit/injection/main.cc).

Another important change has been made to models that consider switchable primary variables. For these models, the `PrimaryVariables` object additionally stores information on the present phases, which means that when defining primary variables, e.g. for Dirichlet boundary conditions or initial conditions, you have to set the phase presence along with the other variables by passing it to the `setState()` function in the `PrimaryVariables`. This made the interfaces `initialPhasePresence()` and `initialPhasePresenceAtPos()` in the __problem__ class obsolete and they have been removed.

__Initial conditions__:

The signature of the interface for defining initial conditions, i.e.

```cpp
void initial(PrimaryVariables &values,
             const Element &element,
             const FVElementGeometry &fvGeometry,
             const int scvIdx) const
```

was changed to

```cpp
template<class Entity>
PrimaryVariables initial(const Entity& entity) const
```

Here, `Entity` can be either a grid vertex (for the box scheme) or a grid element (cc schemes). Note that the function now returns an object of `PrimaryVariables` instead of receiving a reference.

__Boundary conditions__:

The signatures of the interfaces for defining boundary conditions, i.e.

```cpp
void boundaryTypes(BoundaryTypes &values, const Intersection &is) const // cc scheme
void boundaryTypes(BoundaryTypes &values, const Vertex &v) const // box scheme
```

were changed to

```cpp
BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const // cc schemes
BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const // box scheme
```

As for the initial conditions, the functions now return an object of `BoundaryTypes` instead of receiving a reference. This also holds for the specification of boundary condition values, i.e. the interfaces for Neumann or Dirichlet boundary conditions, which have been changed to the following signatures accordingly:

```cpp
PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const // cc schemes
PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const // box scheme

NumEqVector neumann(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars,
                    const SubControlVolumeFace& scvf) const

```

Note that the return type of the `neumann()` function is now of type `NumEqVector`, which might differ from the `PrimaryVariables` type in the fact that it doesn't hold any additional information such as e.g. the `state` in switchable primary variables. It should also be mentioned that the corresponding functions with the suffix `atPos`, taking only a position as argument, still exist and can be used in the same way as in DuMuX 2.12.

Please also note that it is no longer possible to set mixed boundary conditions for cell-centered schemes as this cannot be implemented in a general way in a satisfactory manner. Also, `outflow` boundary conditions are no longer supported, for the same reasons. However, you can still achieve both of the above mentioned boundary conditions manually within the body of the `neumann()` function (see e.g. [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/test/porousmediumflow/1pnc/implicit/1p2c/isothermal/main.cc)). Values returned by the `neumann()` function can depend on primary variables in any fashion.

## Spatial params class
The base classes for __spatial parameters__ in DuMuX 2.12, i.e. `ImplicitSpatialParamsOneP` and `ImplicitSpatialParams`, were moved to the classes `FVSpatialParamsOneP` and `FVSpatialParams`. The classes explicitly take the template parameters `FVGridGeometry` and `Scalar` instead of `TypeTag` in 2.12. The interfaces for returning a spatial parameter (e.g. _porosity_, _permeability_ or _material law parameters_) have been changed from (here shown for the _porosity_)

```cpp
void porosity(const Element &element,
              const FVElementGeometry &fvGeometry,
              int scvIdx) const
```

to

```cpp
template<class ElementSolution>
Scalar porosity(const Element& element,
                const SubControlVolume& scv,
                const ElementSolution& elemSol) const
```

The additional parameter element solution contains the primary variables on all degrees of freedom that are embedded in this element and can be used to realize solution-dependent parameters. Note that the interface name for the parameter _permeability_ has been changed from `intrinsicPermeability()` to `permeability()`.

The wettability of the porous medium can now vary in space/time and/or depending on the current solution. In order to realize this, the `FVSpatialParams` are now equipped with a function that returns the index of the wetting phase within a given `FluidSystem`:

```cpp
template<class FluidSystem, class ElementSolution>
int wettingPhase(const Element& element,
                 const SubControlVolume& scv,
                 const ElementSolution& elemSol) const

template<class FluidSystem>
int wettingPhaseAtPos(const GlobalPosition& globalPos) const
```

The `MaterialLaw` is no longer a dumux _property_, which is why you now have to provide the public aliases `MaterialLaw` and `MaterialLawParams` in your `SpatialParameters` implementation. Furthermore, you must export the type you are using for the _permeability_ (e.g. a scalar or tensor):

```cpp
using MaterialLaw = RegularizedVanGenuchten<Scalar>;
using MaterialLawParams = typename MaterialLaw::Params;
using PermeabilityType = Scalar;
```