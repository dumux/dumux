# Grid

Grids are provided by DUNE or external implementations, but all have to share the DUNE grid interface.
The grid is an abstract concept consisting of several entities: vertices, edges, faces, cells, elements,... Furthermore, additional properties characterize the capabilities and specializations of a grid:
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
* [FoamGrid](https://gitlab.dune-project.org/extensions/dune-foamgrid): 1D/2D grids in physical space of arbitrary dimension, no manifold structures expected. May be used for simulating foams, discrete fracture networks, network flow problems
* [SPGrid](https://gitlab.dune-project.org/extensions/dune-spgrid) (Sparse Paged Grid): structured, parallel
* [MMESHGrid](https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh): can handle moving, physical interfaces
* [SubGrid](https://gitlab.dune-project.org/extensions/dune-subgrid): allows to mark subset of a grid's elements. The marked elements can then be treated in a separate hierarchy.
* OPMGrid (Open Porous Media Grid): supports corner-point format
* [Many more](https://www.dune-project.org/groups/grid/)

Within Dumux, the grid type is set within the `properties.hh` file. Standard grids like YASPGrid, AlbertaGrid and OneDGrid are included in the dune-grid core module. Using other grid types requires downloading the respective DUNE module. The [base implementation](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/grid.hh) of the grid data member is part of the dune-grid module.


### Key functionalites
* maxLevel():
    - Return maximum level defined in this grid. The base grid has level zero, each refinement increases the level by one.
* size(...):
    - Return number of grid or leaf entities of a given `codim` on a given `level`. Via `codim`, which is the complement to the dimension `dim`, geometrical entities can be desribed. If we take a 3D element as an example, then the element itself is characterized by a dimension `dim` of 3 but its `codim` is 0. A face of such element would have `dim` 2 and but its `codim` is 1. An edge of this element would have `dim=1` and `codim=2`. Basically, the `codim` of a point entity is always 0 while a grid element has full dimesion `dim`.
* levelGridView(level):
    - View for a grid `level`. Contains information on all the elements belonging to one certain hierarchy `level` of a grid.
* leafGridView():
    - View for the leaf grid. Contains information on all leaf elements, this is the highest refinement level of a grid.
* globalRefine(refCount):
    - Refine the grid `refCount` times using the default refinement rule.
* mark(refCount, entity):
    - Marks an entity to be refined/coarsened in a subsequent adapt.
* preAdapt():
    - To be called after entities have been marked and before `adapt()` is called. Can include custom operations which should be executed before adapting such as processing or storing data.
* adapt():
    - Refine all positive marked leaf entities, coarsen all negative marked entities if possible.
* postAdapt():
    - To be called after grid has been adapted and information left over by the adaptation has been processed.
* comm():
    - Return the communication object which manages parallelization in case of a parallel grid.
* loadBalance():
    - Re-balances the load each process has to handle for a parallel grid.

### Overview

@mermaid{grid}