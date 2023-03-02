## grid

Grids are provided by DUNE or external implementations, but all have to share the DUNE grid interface.
The grid is an abstract concept and consists of several entities: vertices, edges, faces, cells, elements,... Furthermore, additional properties characterize the capabilities and specializations of a grid:
* dimension and co-dimension
* element types
* conforming or nonconforming
* hierarchical
* local refinement, refinement rules
* parallel data distribution and communication, dynamic load balancing.

Depending on the needs, one has to choose the proper grid. Some grids and what characterizes them are:
* [YASPGrid](https://gitlab.dune-project.org/core/dune-grid) (Yet Another Structured Parallel Grid): structured, parallel, arbitrary overlap
* [AlbertaGrid](https://gitlab.dune-project.org/core/dune-grid): provides access to ALBERTA finite element toolbox via DUNE interface
* [OneDGrid](https://gitlab.dune-project.org/core/dune-grid) (One Dimensional Grid): 1D, adaptive
* [ALUGrid](https://gitlab.dune-project.org/extensions/dune-alugrid) (adaptive, loadbalancing, unstructured Grid): adaptive, loadbalancing, unstructured
* [UGGrid](https://gitlab.dune-project.org/staging/dune-uggrid) (Unstructured Grid): 2D/3D, unstructured grid, can be refined
* [FoamGrid](https://gitlab.dune-project.org/extensions/dune-foamgrid): 1D/2D grids in physical space of arbitrary dimension, no manifold structures expected. Best for simulating foams, discrete fracture networks, network flow problems
* [SPGrid](https://gitlab.dune-project.org/extensions/dune-spgrid) (Sparse Paged Grid): structured, parallel
* [MMESHGrid](https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh): can handle moving, physical interfaces
* [SubGrid](https://gitlab.dune-project.org/extensions/dune-subgrid): allows to mark subset of a grid's elements and these elements can be treated in a hierarchy of their own
* OPMGrid (Open Porous Media Grid): supports corner-point format
* [Many more](https://www.dune-project.org/groups/grid/)

Within Dumux, the grid type is set within the `properties.hh` file. Standard grids like YASPGrid, AlbertaGrid and OneDGrid are included in the dune-grid core module. Using other grid types requires downloading the respective DUNE module. The [base implementation](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/grid.hh) of the grid data member is part of the dune-grid module.


### Key functionalites
1. maxLevel()
2. size()
3. levelGridView()
4. leafGridView()
5. globalRefine()
6. mark()
7. preAdapt()
8. postAdapt()
9. comm()
10. loadBalance()

### Overview
```{mermaid}
flowchart LR
    A(Grid) -->|"template parameter for declaration"| B(gridManager)
    B -->|"init()"| C(grid)
    click A "./grid.html"
    click B "./gridmanager.html"
    click C "./grid.html"
```
