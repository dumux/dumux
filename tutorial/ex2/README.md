# Exercise #2 (DuMuX course)

## Problem set-up

The problem setup is identical to the previous [_exercise 1_](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/6c2ca5dd000599f793f89fe5e3fc69ef9c9d8b73/tutorial/ex1/README.md) with a lower injection rate of 1e-6 kg/(m*s) so that diffusion plays a more dominant role in the transport process.

## Preparing the exercise

* Navigate to the directory `dumux/tutorial/dumux-course`

_Exercise 2_ deals with a two-phase compositional problem (__2p2c__). Goal is to learn how to use compile and runtime parameters and the _DuMuX property system_.

### 1. Getting familiar with the code

Locate all the files you will need for this exercise
* The __main file__: `exercise2.cc`
* The __problem file__: `injection2p2cproblem.hh`
* The __spatial parameters file__: `injection2pspatialparams.hh`
* The __input file__: `exercise2.input`
* Two header files containing:
  * a custom __local residual__ in: `mylocalresidual.hh`
  * a custom __material law__ in:: `mymateriallaw.hh`

### 2. Compiling and running the program

* Change to the build-directory

```bash
cd build-cmake/tutorial/exercise2
```

* Compile the executable `exercise2`

```bash
make exercise2
```

* Execute the two problems and inspect the result

```bash
./exercise2
```
Note: Because the input file has the same name as the
executable, DuMuX will find it automatically.

If gnuplot is installed on your system, you should see a plot of the capillary pressure - saturation relationship.

### 3. Implement and use a different material law

DuMuX uses the term _material law_ to describe the law used to compute
* pc-Sw relations
* kr-Sw relations
* their inverse relations

The file `mymateriallaw.hh` contains a custom implementation of such a material law.

* Implement the method `Scalar pc(Scalar sw)` by implementing your own capillary pressure relationship, e.g. pc(Sw) = 1e5*(1-Sw).

The type (i.e. C++ type) of the material law is set in the file `injection2pspatialparams.hh` by using the DuMuX property system

```c++
SET_PROP(InjectionSpatialParams, MaterialLaw)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
};
```

* Make DuMuX use your own material law by including the header `mymateriallaw.hh` and changing the alias `type`. This will make sure that your material law is used everywhere else in the code.

Verify your changes by recompiling and running the program. You should see a plot of your new function.

### 4. Implement your own local residual

Most types in DuMuX are properties that can be changed just like the material law. In the following task we implement our own 2p2c local residual, i.e. the class that computes the element residual  in every Newton step. The file `mylocalresidual.hh` contains a copy of the original local residual class used for the 2p2c model renamed to `template<class TypeTag> class MyTwoPTwoCLocalResidual`.

* Make DuMuX use this new local residual by setting the corresponding property in the `Property` namespace in the file `injection2p2cproblem.hh`

```c++
// note that every property struct knows about TypeTag
SET_PROP(Injection2p2cProblem, LocalResidual)
{
    using type = MyTwoPTwoCLocalResidualocal<TypeTag>;
};

// or using the convenience macro
SET_TYPE_PROP(Injection2p2cProblem, LocalResidual,
              MyTwoPTwoCLocalResidualocal<TypeTag>);
```

* Implement an output to the terminal in the constructor of `MyTwoPTwoCLocalResidual` e.g.

```c++
MyTwoPTwoCLocalResidual()
{
    std::cout << "Using MyTwoPTwoCLocalResidual." << std::endl;
}
```

* Verify you are using the new class by compiling and running the new program and inspecting the terminal output.

You want to make the new local residual special by adding a switch enabling / disabling diffusion. We will achieve this with a DuMuX parameter, a parameter read from the input file that defaults to a property value if the input file doesn't contain the parameter.

* Create a new `TypeTag` node, a new `PropertyTag`, and set a default in the `mylocalresidual.hh` file by adding

```c++
namespace Dumux {

namespace Properties
{
    NEW_TYPE_TAG(MyLocalResidualParams); // creates a new TypeTag node
    NEW_PROP_TAG(ProblemEnableDiffusion); // creates a new property
    SET_BOOL_PROP(MyLocalResidualParams,
                  ProblemEnableDiffusion, true); // set a default value
}
...
```

* Modify the `computeFlux` method to only call the `diffusiveFlux` method if diffusion is enabled. You can get the new parameter by adding the lines
```c++
// ... in the constructor of MyTwoPTwoCLocalResidual
    enableDiffusion_ = GET_PARAM_FROM_GROUP(TypeTag, bool,
                            Problem, EnableDiffusion);

// ... in the private member section of MyTwoPTwoCLocalResidual
private:
    bool enableDiffusion_;
```
