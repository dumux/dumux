# Exercise #1 (DuMuX course)

## Problem set-up

N2 is injected in an aquifer previously saturated with water with an injection rate of  $`0.001~kg/m*s`$.
The aquifer is situated 2700 m below see level and the domain size is 60 m x 40 m. It consists of two layers, a moderately permeable one ($`\Omega 1`$) and a lower permeable one ($`\Omega 2`$).

<img scr="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/raw/feature/dumux-course-exercise1/tutorial/extradoc/exercise1_setup.png" width="200"/>

## Preparing the exercise

* Navigate to the directory `dumux/tutorial/dumux-course`

_Exercise 1_ deals with two problems: a two-phase immiscible problem (__2p__) and a two-phase compositional problem (__2p2c__). They both set up the same scenario with the difference that the 2p2c assumes a miscible fluid state for the two fluids (water and gaseous N2) and the 2p model assumes an immiscible fluid state.

### 1. Getting familiar with the code

Locate all the files you will need for this exercise
* The __main file__ for the __2p__ problem: `exercise1_2p.cc`
* The __problem file__ for the __2p__ problem: `injection2pproblem.hh`
* The __main file__ for the __2p2c__ problem: `exercise1_2p2c.cc`
* The problem file for the __2p2c__ problem: `injection2p2cproblem.hh`
* The shared __spatial parameters file__: `injection2pspatialparams.hh`
* The shared __input file__: `exercise1.input`


### 2. Compiling and running an executable

* Change to the build-directory

```bash
cd build-cmake/tutorial/exercise1
```

* Compile both executables `exercise1_2p` and `exercise1_2p2c`

```bash
make exercise1_2p exercise1_2p2c
```

* Execute the two problems and inspect the result

```bash
./exercise1_2p exercise1.input
./exercise1_2p2c exercise1.input
```

### 3. Changing input parameters

In the input file `exercise1.input` you can find the following section

```ini
[SpatialParams]
EntryPressureFine = 1e4
EntryPressureCoarse = 1e4
```

* Change the values for the fine entry pressure in the input file to a lower value and compare the results with the previous solution. You don't need to recompile the executable.

### 4. Runtime parameters

The injection rate is currently hard-coded in `injection2pproblem.hh`.

```c++
code snippet
```

We want to be able to set it at runtime. To this end,
* use the following DuMuX macro to read a runtime parameter from the input file

```c++
// read the injection rate from the input file at run time
const auto injectionRate = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, <TYPE>, <GROUPNAME>, <PARAMNAME>);
```

* Replace
`<TYPE>`,`<GROUPNAME>`,`<PARAMNAME>` by what is appropriate for your case:
  * `<TYPE>` is the type of the parameter to read
  * `<GROUPNAME>` is the group in the input file
  * `<PARAMNAME>` is the name of the parameter in the input file

Note that due to the way the macro works, the names are specified as plain text without "quotation marks".

* Check the influence of that parameter on the simulation result by rerunning the simulation with different injection rates.

Again, you don't need to recompile the program.

### 5. Setting up a new executable (for a non-isothermal simulation)

* Set up a new cc file called `exercise1_2pni.cc` by copying and renaming `exercise1_2p.cc`

```bash
cp exercise1_2p.cc exercise1_2pni.cc
```

* In  `exercise1_2pni.cc`, include the header `injection2pniproblem.hh` instead of the isothermal problem file `injection2pproblem.hh`.
* Add a new executable in `CMakeLists.txt` by adding the lines

```cmake
dune_add_test(NAME injection2pniproblem
              SOURCES injection2pniproblem.cc)
```

* Test that everything compiles without error

```bash
make # should rerun cmake
make injection2pniproblem # builds new executable
```

### 6. Setting up a non-isothermal __2pni__ test problem

* Open the file `injection2pniproblem.hh`. It is a copy of the `injection2pproblem.hh` with some useful comments on how to implement a non-isothermal model. Look for comments containing

```c++
// TODO: dumux-course-task
```

* The following set-up should be realized:

  __Boundary conditions:__ Dirichlet conditions for the temperature with a temperature gradient of 0.03 K/m and a starting temperature of 283 K.

  __Initial conditions:__ The same temperature gradient as in the boundary conditions with an additional lens (with position: 20 < x < 30, 5 < y < 35), which has an initial temperature of 380 K.

<img scr="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/raw/feature/dumux-course-exercise1/tutorial/extradoc/exercise1_nonisothermal.png" width="200"/>

The non-isothermal model requires additional spatial parameters like the thermal conductivity. They are already implemented in `injection2pspatialparams.hh`, you just need to _uncomment_ them.
