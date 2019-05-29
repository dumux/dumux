This tutorial was copied from tests/porousmediumflow/2p/adaptive and restricted to the cell-centered finite volume TPFA discretization scheme.
You need [ALUGrid][0] in order to compile and run it.

# Two-phase flow with infiltration and adaptive grid

## Problem set-up
This example models a soil contamination problem where DNAPL infiltrates a fully water saturated medium.

The two-dimensional model domain is 6m x 4m and contains a lens with a lower permeability for 1m ≤ x ≤ 3m and 2m ≤ y ≤ 4m.
A linear pressure gradient is given as a Dirichlet boundary condition at the left and the right boundary.
Neumann boundary conditions are set at the upper and lower boundary.
DNAPL enters the model domain at the upper boundary for 1.75m ≤ x ≤ 2m with a rate of 0.04kg/ms.
In addition, the DNAPL is injected at a point source at x = 0.502 and y = 3.02 with a rate of 0.1kg/s.


## Infiltration (point source)
### problem.hh
The point sources are specified in the problem.hh file by the `addPointSources` method.
You can add an arbitrary number of point sources to the vector of point sources.
To instantiate a point source, the position and the infiltration values are needed.
For the definition of the `PointSource` class see dumux/common/pointsource.hh
```C++
void addPointSources(std::vector<PointSource>& pointSources) const
{
    // inject 0.1 kg/s of non-wetting phase at position (0.502, 3.02);
    pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
}
```

### main.cc
The `computePointSourceMap` method is called from the main.cc file to compute the point sources. It must be called during the initialisation (l. 97) and after each refinement of the mesh (ll. 135, 151 and 205).
```C++
problem->computePointSourceMap();
```
The `computePointSourceMap` method is inherited from the fvproblem and therefore specified in the dumux/common/fvproblem.hh.
It calls the `addPointSources` method specified in the problem.hh file.

## Adaptive grid
### main.cc
The following headers need to be included in the main file:

```C++
#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/adaptive/initializationindicator.hh>
#include <dumux/porousmediumflow/2p/griddatatransfer.hh>
#include <dumux/porousmediumflow/2p/gridadaptindicator.hh>
```

The grid adaptation is prepared during the initialization by the following steps:
1. **Instantiate indicator (l. 114):**
The indicator is saturation-dependent and defined in the file dumux/porousmediumflow/2p/gridadaptindicator.hh.
It allows to set the minimum and maximum allowed refinement levels via the input parameters
`Adaptive.MinLevel` and `Adaptive.MaxLevel`.
2. **Instantiate data transfer (l. 154):**
The data transfer performs the transfer of data on a grid from before to after adaptation and is defined in the
file dumux/porousmediumflow/2p/griddatatransfer.hh.
Its main functions are to store and reconstruct the primary variables.
3. **Set the indicator for the initial refinement around sources/BCs (l. 118):**
We use the GridAdaptInitializationIndicator defined in dumux/adaptive/initializationindicator.hh.
4. **Refine up to the maximum level (ll. 121-137):**
For every level, the indicator used for the refinement/coarsening is calculated.
If any grid cells have to be adapted, the gridvariables and the pointsourcemap are updated.
5. **Do refinement for the initial conditions using the indicator (ll. 140-144):**
Depending on the initial conditions, another grid adaptation might be necessary.
The gridadaptindicator uses the input parameters `Adaptive.RefineTolerance` and `Adaptive.CoarsenTolerance` for this step.
For further details on the indicator calculations see dumux/porousmediumflow/2p/gridadaptindicator.hh, ll. 116.
Afterwards, the marked elements are adapted.
6. **Update grid data after adaption (ll. 147-152):**
In case of a grid adaptation, the gridvariables and the pointsourcemap are updated.


During the time loop, the refinement indicator is computed (l. 191) and the respective elements to be refined are marked (ll. 194-196).

In case of grid adaptation, the following updates are necessary (ll. 201-205):
```C++
 xOld = x; //!< Overwrite the old solution with the new (resized & interpolated) one
 ssembler->setJacobianPattern(); //!< Tell the assembler to resize the matrix and set pattern
 assembler->setResidualSize(); //!< Tell the assembler to resize the residual
 gridVariables->updateAfterGridAdaption(x); //!< Initialize the secondary variables to the new (and "new old") solution
 problem->computePointSourceMap(); //!< Update the point source map
```

### problem.hh
A non-conforming grid such as ALUGrid has to be chosen in the problem file:
```C++
  //! Use non-conforming refinement
  template<class TypeTag>
  struct Grid<TypeTag, TTag::TwoPAdaptivePointSource> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
```


### params.input
The following parameters in the `[Adaptive]` parameter group determine the grid adaptation behavior:
* `RefineAtDirichletBC`: If to refine at Dirichlet boundaries
* `RefineAtFluxBC`: If to refine at Neumann/Robin boundaries
* `MinLevel`: Minimum allowed refinement level, used by the indicators
* `MaxLevel`: Maximum allowed refinement level, used by the indicators
* `CoarsenTolerance`
* `RefineTolerance`

## Solution
![](test_2p_pointsource_adaptive.png)


[0]: https://gitlab.dune-project.org/extensions/dune-alugrid
