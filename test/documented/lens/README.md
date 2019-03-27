# Richards lense setup using tensor product coordinates

## Problem set-up

 The domain is box shaped. Left and right boundaries are Dirichlet boundaries with fixed water pressure (fixed Saturation `S_w = 0`),  bottom boundary is closed (Neumann 0 boundary), the top boundary (Neumann 0 boundary) is also closed except for infiltration section, where water is infiltrating into an initially unsaturated porous medium. This problem is very similar to the LensProblem which uses the TwoPBoxModel, with the main difference being that the domain is initially fully saturated by gas instead of water and water instead of a DNAPL infiltrates from the top.

## Tensor product coordinates

### What is it used for
* enables to grade the mesh distribution.

### What to change code
`problem.hh`

**original code:**
```c++
struct Grid<TypeTag, TTag::RichardsLens> { using type = Dune::YaspGrid<2>; };
```

**code including the tensor product coordinate feature:**

```c++
struct Grid<TypeTag, TTag::RichardsLens> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };
```

`params.input`

**original code:**
```c++
[Grid]
UpperRight = 6 4
Cells = 24 16
```
**code including the tensor product coordinate feature:**
```c++
[Grid]
Cells0 = 24
Cells1 = 30
Positions0 = 0.0 6
Positions1 = 0.0 4
Grading0 = 1.0
Grading1 = -1.2
```
