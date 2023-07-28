# Grid

Grids represent the spatial discretization of a physical domain in an abstract way by approximating the domain via several entities: vertices, edges, faces, cells, elements... These entities can be accesses via the *grid view* object in Dumux, see @ref gridview.
While grids can be provided by DUNE or external implementations, all grid types have to share a common DUNE grid interface.

Additional properties characterize the capabilities and specializations of a grid:
* dimension and co-dimension
* element types
* conforming or nonconforming
* hierarchical
* local refinement, refinement rules
* parallel data distribution and communication, dynamic load balancing.

Depending on the needs, one has to choose the proper grid. The core module `dune-grid` already offers a few options, some of these are:
* [YASPGrid](https://gitlab.dune-project.org/core/dune-grid) (Yet Another Structured Parallel Grid): structured, parallel, arbitrary overlap.
* [AlbertaGrid](https://gitlab.dune-project.org/core/dune-grid): provides access to ALBERTA finite element toolbox via DUNE interface.
* [OneDGrid](https://gitlab.dune-project.org/core/dune-grid) (One Dimensional Grid): 1D, adaptive.
* There are other external grid modules which specialize in certain aspects of grid capabilites. For an overview, see @ref external-libraries.


Within Dumux, the grid type is set within the `properties.hh` file. Standard grids like YASPGrid, AlbertaGrid and OneDGrid are included in the dune-grid core module. Using other grid types requires downloading the respective DUNE module. The [base implementation](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/grid.hh) of the grid data member is part of the dune-grid module. Detailed information about DUNE and its modules can be found in the book [DUNE - The Distributed and Unified Numerics Environment](https://link.springer.com/book/10.1007/978-3-030-59702-3#bibliographic-information) by Oliver Sanders.

### Example
```cpp
...
using Grid = GetPropType<TypeTag, Properties::Grid>; //defined in properties as e.g. 'using Grid = Dune::YaspGrid<2>;'
GridManager<Grid> gridManager;
gridManager.init();
...
```
The functionalites of the grid itself are usually not used directly. The `Grid` type is rather used as a template argument for creating a @ref gridmanager "gridManager" object.


### Key functionalites
* See [base implementation](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/grid.hh).

### Overview

@mermaid{grid}
