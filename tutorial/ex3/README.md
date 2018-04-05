# Exercise #3 (DuMuX course)

The aim of this exercise is to get familiar with the _DuMuX_ way of implementing new components (fluids) and fluid systems (mixtures). In the scope of this exercise, a new fictitious component is implemented (exercise _3a_) as well as its mixture with water (exercise _3b_).

## Problem set-up

The domain has a size of 60 x 60 m and contains two low-permeable lenses. Initially, the domain is fully water saturated and the fictitious component is injected through the middle portion of the upper boundary by means of a Neumann boundary condition. The remaining parts of the upper and the entire lower boundary are Neumann no-flow while on the two lateral sides Dirichlet boundary conditions are applied (hydrostatic conditions for the pressure and zero saturation).

![](../extradoc/exercise3_setup.png)


## Preparing the exercise

* Navigate to the directory `dumux/tutorial/ex3`

### 1. Getting familiar with the code

Locate all the files you will need for this exercise
* The shared __main file__ : `exercise3.cc`
* The __input file__ for part a: `exercise3_a.input`
* The __problem file__ for part a: `2pproblem.hh`
* The __input file__ for part b: `exercise3_b.input`
* The __problem file__ for part b: `2p2cproblem.hh`
* The __spatial parameters file__: `spatialparams.hh`

Furthermore you will find the following folders:
* `binarycoefficients`: Stores headers containing data/methods on binary mixtures
* `components`: Stores headers containing data/methods on pure components
* `fluidsystems`: Stores headers containing data/methods on mixtures of pure components. Uses methods from `binarycoefficients`.

To see more components, fluidsystems and binarycoefficients implementations, have a look at the folder `dumux/material`.

### 2. Implement a new component

In the following the basic steps required to set the desired fluid system are outlined. Here, this is done in the __problem file__, i.e. for this part of the exercise the code shown below is taken from the `2pproblem.hh` file.

In this part of the exercise we will consider a system consisting of two immiscible phases. Therefore, the _TypeTag_ for this problem derives from the _BoxTwoP_ _TypeTag_ (for a cell-centered scheme, you would derive from _CCTwoP_). In order to be able to derive from this _TypeTag_, The declaration of the _BoxTwoP_ _TypeTag_ is found in the file `dumux/porousmediumflow/2p/implicit/model.hh`, which has been included in line 28:

```c++
// The numerical model
#include <dumux/porousmediumflow/2p/implicit/model.hh>
```

 Additionally we derive from the _ExerciseThreeSpatialParams_ _TypeTag_ in order to set all the types related to the spatial parameters also for our new _TypeTag_ _ExerciseThreeProblem_. In this case, only one property is set for _ExerciseThreeSpatialParams_ (see lines 43 to 60 in `spatialparams.hh`). Alternatively, we could have defined this also here in the problem without defining a new type tag for the spatial parameters. However, in case you want to define several properties related to the spatial parameters it is good practice to define a separate _TypeTag_ and derive from this, as it is done here.

```c++
NEW_TYPE_TAG(ExerciseThreeTypeTag, INHERITS_FROM(BoxTwoP, ExerciseThreeSpatialParams));
```

As wetting phase we want to use water and we want to precompute tables on which the properties are then interpolated in order to save computational time. Thus, in a first step we have to include the following headers:

```c++
// The water component
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
```
The non-wetting phase will be our new component, where we want to implement an incompressible and a compressible variant. The respective headers are prepared, but still incomplete. The compressible variant is still commented so that compilation does not fail when finishing the incompressible variant.

```c++
// The components that will be created in this exercise
#include "components/myincompressiblecomponent.hh"
// #include "components/mycompressiblecomponent.hh"
```
As mentioned above, we want to simulate two non-mixing components. The respective fluid system is found in:

```c++
// The two-phase immiscible fluid system
#include <dumux/material/fluidsystems/2pimmiscible.hh>
```

This fluid system expects __phases__ as input and so far we have only included the components, which contain data on the pure component for all physical states. Thus, we need to include

```c++
// We will only have liquid phases Here
#include <dumux/material/fluidsystems/1pliquid.hh>
```

which creates a _liquid phase_ from a given component. Finally, using all of the included classes we set the the fluid system property by choosing a liquid phase consisting of the incompressible fictitious component as non-wetting phase and tabulated water as the wetting phase in the immiscible fluid system:


```c++
// we use the immiscible fluid system here
SET_PROP(ExerciseThreeTypeTag, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Components::TabulatedComponent<Components::H2O<Scalar>> TabulatedH2O;
    typedef typename FluidSystems::OnePLiquid<Scalar, TabulatedH2O> WettingPhase;
    /*!
     * Uncomment first line and comment second line for using the incompressible component
     * Uncomment second line and comment first line for using the compressible component
     */
    typedef typename FluidSystems::OnePLiquid<Scalar, MyIncompressibleComponent<Scalar> > NonWettingPhase;
    // typedef typename FluidSystems::OnePLiquid<Scalar, MyCompressibleComponent<Scalar> > NonWettingPhase;

public:
    typedef typename FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonWettingPhase> type;
};
```

### 2.1. Incompressible component

Open the file `myincompressiblecomponent.hh`. You can see in line 40 that a component should always derive from the _Base_ class (see `dumux/material/components/base.hh`), which defines the interface of a _DuMuX_ component with all possibly required functions to be overloaded by the actual implementation. Additionally it is required for liquids to derive from the _Liquid_ class (see `dumux/material/components/liquid.hh`), for gases to derive from the _Gas_ class (see `dumux/material/components/gas.hh`) and for solids to derive from the _Solid_ class (see `dumux/material/components/solid.hh`).

```c++
/*!
 * \ingroup Components
 * \brief A ficitious component to be implemented in exercise 3.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class MyIncompressibleComponent : public Component<Scalar, MyIncompressibleComponent<Scalar> >
```

__Task__:

Implement an incompressible component into the file `myincompressiblecomponent.hh`, which has the following specifications:

| Parameter | unit | value |
| -----| --------| -------- |
| $`M`$ | $`Kg/mol`$   | $`131.39 \cdot 10^{-3}`$ |
| $`\rho_{liquid}`$ | $`Kg/m^3`$   | $`1460`$   |
| $`\mu_{liquid}`$ | $`Pa \cdot s`$   | $`5.7 \cdot 10^{-4}`$   |

In order to do so, have a look at the files `dumux/material/components/base.hh` and `dumux/material/components/liquid.hh` to see how the interfaces are defined and overload them accordingly.

In order to execute the program, change to the build-directory

```bash
cd build-cmake/tutorial/ex3
```

Uncomment the line for the corresponding executable in the `CMakeLists.txt` file:

```cmake
dune_add_test(NAME exercise3_a
              SOURCES exercise3.cc
              COMPILE_DEFINITIONS TYPETAG=ExerciseThreeBoxTwoPTypeTag
              CMD_ARGS exercise3_a.input)
```

Now you can compile and execute the program by typing

```bash
make
make exercise3_a
./exercise3_a exercise3_a.input
```

The saturation distribution at the final simulation time should look like this:

![](../extradoc/exercise3_a_solution.png)

### 2.2. Compressible component

We now want to implement a pressure-dependent density for our component. Open the file `mycompressiblecomponent.hh` and copy in the functions you implemented for the incompressible variant. Now substitute the method that returns the density by the following expression:

$`\displaystyle \rho_{MyComp} = \rho_{min} + \frac{ \rho_{max} - \rho_{min} }{ 1 + \rho_{min}*e^{-1.0*k*(\rho_{max} - \rho_{min})*p} } `$

where $`p`$ is the pressure and $`\rho_{min} = 1440 `$, $`\rho_{max} = 1480 `$ and $`k = 5 \cdot 10^{-7} `$. Also, make sure the header is included in the `2pproblem.hh` file by uncommenting line 42. Furthermore, the new component has to be set as the non-wetting phase in the fluid system, i.e. comment line 81 and uncomment line 82. The non-wetting density distribution at the final simulation time should look like this:

![](../extradoc/exercise3_a_solution2.png)

### 3. Implement a new fluid system

The problem file for this part of the exercise is `2p2cproblem.hh`. We now want to implement a new fluid system consisting of two liquid phases, which are water and the previously implemented compressible component. We will consider compositional effects, which is why we now have to derive our _TypeTag_ from the _BoxTwoPTwoC_ _TypeTag_:

```c++
// The numerical model
#include <dumux/porousmediumflow/2p2c/model.hh>
```

```c++
// Create a new type tag for the problem
NEW_TYPE_TAG(ExerciseThreeTypeTag, INHERITS_FROM(BoxTwoPTwoC, ExerciseThreeSpatialParams));
```

The new fluid system is to be implemented in the file `fluidsystems/h2omycompressiblecomponent.hh`. This is already included in the problem and the fluid system property is set accordingly.

```c++
// The fluid system that is created in this exercise
#include "fluidsystems/h2omycompressiblecomponent.hh"
```

```c++
// The fluid system property
SET_PROP(ExerciseThreeTypeTag, FluidSystem)
{
private:
   typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
   typedef FluidSystems::H2OMyCompressibleComponent<TypeTag, Scalar> type;
};
```

In the `fluidsystems/h2omycompressiblecomponent.hh` file, your implemented component and the binary coefficient files are already included.

```c++
// the ficitious component that was created in exercise 3a
#include <tutorial/ex3/components/mycompressiblecomponent.hh>

// the binary coefficients corresponding to this fluid system
#include <tutorial/ex3/binarycoefficients/h2omycompressiblecomponent.hh>
```

__Task__:

Under the assumption that one molecule of _MyCompressibleComponent_ displaces exactly one molecule of water, the water phase density can be expressed as follows:

$` \rho_{w} = \frac{ \rho_{w, pure} }{ M_{H_2O} }*(M_{H_2O}*x_{H_2O} + M_{MyComponent}*x_{MyComponent}) `$

Implement this dependency in the _density()_ method in the fluid system. In order to compile and execute the program, uncomment the line for the corresponding executable in the `CMakeLists.txt` file

```cmake
dune_add_test(NAME exercise3_b
              SOURCES exercise3.cc
              COMPILE_DEFINITIONS TYPETAG=ExerciseThreeBoxTwoPTwoCTypeTag
              CMD_ARGS exercise3_b.input)
```

Then, change to the build-directory

```bash
cd build-cmake/tutorial/ex3
```

and type

```bash
make
make exercise3_b
./exercise3_b exercise3_b.input
```

You will observe an error message and an abortion of the program. This is due to the fact that in order for the constraint solver and other mechanisms in the two-phase two-component model to work, two additional functionalities in the component have to be implemented. The model has to know whether or not the liquid pure component is compressible and it needs the vapour pressure. As in the previous exercise, check the `dumux/material/components/base.hh` and `dumux/material/components/liquid.hh` files for these two functions and implement them into `mycompressiblecomponent.hh`. For the vapour pressure, use a value of $`3900`$  Pa.
