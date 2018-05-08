# Exercise #1 (DuMuX course)
<br>
## Problem set-up

N$`_2`$ is injected in an aquifer previously saturated with water with an injection rate of 0.001 kg/(s*m$`^2`$).
The aquifer is situated 2700 m below sea level and the domain size is 60 m x 40 m. It consists of two layers, a moderately permeable one ($`\Omega_1`$) and a lower permeable one ($`\Omega_2`$).

<img src="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/raw/master/tutorial/extradoc/exercise1_setup.png" width="1000">

## Preparing the exercise

* Navigate to the directory `dumux/tutorial/ex1`

_Exercise 1_ deals with two problems: a two-phase immiscible problem (__2p__) and a two-phase compositional problem (__2p2c__). They both set up the same scenario with the difference that the 2p2c assumes a miscible fluid state for the two fluids (water and gaseous N$`_2`$) and the 2p model assumes an immiscible fluid state.

<br><br>
### Task 1: Getting familiar with the code
<hr>

Locate all the files you will need for this exercise
* The __main file__ for the __2p__ problem : `exercise1_2p.cc`
* The __main file__ for the __2p2c__ problem : `exercise1_2p2c.cc`
* The __problem file__ for the __2p__ problem: `injection2pproblem.hh`
* The __problem file__ for the __2p2c__ problem: `injection2p2cproblem.hh`
* The shared __spatial parameters file__: `injection2pspatialparams.hh`
* The shared __input file__: `exercise1.input`

<hr><br><br>
### Task 2: Compiling and running an executable
<hr>

* Change to the build-directory

```bash
cd ../../build-cmake/tutorial/ex1
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

* you can look at the results with paraview

```bash
paraview injection-2p2c.pvd
```

<hr><br><br>
### Task 3: Changing input parameters
<hr>

In the input file `exercise1.input` you can find the following section

```ini
[SpatialParams]
EntryPressureAquitard = 4.5e4
EntryPressureAquifer = 1e4
```

* Change the values for the aquitard entry pressure in the input file to a lower value and compare the results with the previous solution. You don't need to recompile the executable.

<hr><br><br>
### Task 4: Runtime parameters
<hr>

The injection rate is currently hard-coded in `injection2p2cproblem.hh` to $`1e-4 kg/(s m^2)`$.

```c++
 // set the Neumann values for the Nitrogen component balance
 // convert from units kg/(s*m^2) to mole/(s*m^2)
values[Indices::contiNEqIdx] = -1e-4/FluidSystem::molarMass(FluidSystem::nCompIdx);
values[Indices::contiWEqIdx] = 0.0;
```

We want to be able to set it at runtime. To this end,
* use the following DuMuX macro to read a runtime parameter from the input file

```c++
// read the injection rate from the input file at run time
totalAreaSpecificInflow_ = getParam<TYPE>("GROUPNAME.PARAMNAME");
```

* Replace
`<TYPE>`,`<GROUPNAME>`,`<PARAMNAME>` by what is appropriate for your case:
  * `<TYPE>` is the type of the parameter to read
  * `<GROUPNAME>` is the group in the input file
  * `<PARAMNAME>` is the name of the parameter in the input file

Note that due to the way the macro works, the names are specified as plain text within the "quotation marks".`<GROUPNAME>` and `<PARAMNAME>` need to be separated by a dot (.).
Follow the instructions given as a

```c++
// TODO: dumux-course-task
```
in the `injection2p2cproblem.hh` file and also remember to also set the parameter totalAreaSpecificInflow in the input file.

* Check the influence of that parameter on the simulation result by rerunning the simulation with different injection rates. Remember to also set the parameter totalAreaSpecificInflow in the input file.

Since you have changed your header file, you have to recompile the program.

<hr><br><br>
### 5. Setting up a new executable (for a non-isothermal simulation)
<hr>

* Copy the main file `exercise1_2p.cc` and rename it to `exercise1_2pni.cc`
* In  `exercise1_2pni.cc`, include the header `injection2pniproblem.hh` instead of `injection2pproblem.hh`.
* In  `exercise1_2pni.cc`, change `Injection2pCCTypeTag` to `Injection2pNICCTypeTag` in the line `using TypeTag = TTAG(Injection2pCCTypeTag);`
* Add a new executable in `CMakeLists.txt` by adding the lines

```cmake
# the two-phase non-isothermal simulation program
dune_add_test(NAME exercise1_2pni
              SOURCES exercise1_2pni.cc
              CMD_ARGS exercise1.input)
```

* Test that everything compiles without error

```bash
make # should rerun cmake
make exercise1_2pni # builds new executable
```

<hr><br><br>
### 6. Setting up a non-isothermal __2pni__ test problem
<hr>

* Open the file `injection2pniproblem.hh`. It is a copy of the `injection2pproblem.hh` with some useful comments on how to implement a non-isothermal model. Look for comments containing

```c++
// TODO: dumux-course-task
```

* The following set-up should be realized:

  __Boundary conditions:__ Dirichlet conditions for the temperature with a temperature gradient of 0.03 K/m and a starting temperature of 283 K.

  __Initial conditions:__ The same temperature gradient as in the boundary conditions with an additional lens (with position: 20 < x < 30, 5 < y < 35), which has an initial temperature of 380 K.

<img src="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/raw/master/tutorial/extradoc/exercise1_nonisothermal.png" width="800">

The non-isothermal model requires additional parameters like the thermal conductivity of the solid component. They are already implemented and set in `exercise1.input`, you just need to _uncomment_ them.
