# GridManager

The responsibility of the grid manager is to create grids, handle data associated with grid entities and manage their distribution across the processes of parallel computations. The grid manager supports different kinds of parameters from which to construct a grid: (1) the path to a grid file, or (2) a set of parameters from which a structured grid can be created, such as the dimensions of the grid and the number of cells to be used for discretization. The latter is available for all supported grid types, while some restrictions may apply for the construction from grid files. As an example, you cannot construct a structured `Dune::YaspGrid` from an unstructured grid file. Grid manager implementations for a wide range of grid types are available, for instance, for `YaspGrid`, `ALUGrid` and `UGGrid`. For a complete overview, have a look at the folder [`dumux/io/grid`](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/io/grid/).

### Example

```cpp
...
using Grid = Dune::YaspGrid<2>;
GridManager<Grid> gridManager; // creates a GridManager specialized for the chosen `Grid` type (here Dune::YaspGrid<2>)
gridManager.init(); // initializes the grid using parameters from the `params.input` file
const Grid& grid = gridManager.grid(); // get a reference to the constructed grid
...
```
Construction from grid files takes precedence over the construction of structured grids. That is, when the `init()` function is called on a grid manager, it checks if the parameter `Grid.File` is defined, in which case it tries to construct the grid from the given file. If this parameter is not defined, it tries to construct a structured grid from the respective parameters. These include the mandatory parameters `Grid.UpperRight` and `Grid.Cells`, which define the upper right corner of the grid and the number of cells used to discretize each dimension. Depending on the grid type, further parameters may be supported, such as `Grid.LowerLeft`, `Grid.Overlap`, `Grid.Positions`, etc. To learn more about the grid type-specific parameters, have a look at the documentation of specific grid manager implementations.

### Key functionalities

* See Dumux::GridManagerBase

### Overview

@mermaid{gridmanager}
