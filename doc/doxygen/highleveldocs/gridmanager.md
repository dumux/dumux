# GridManager

A grid manager handles the grid data from an input (file) and constructs a grid from that information. Only works for supported grid types: YASP, OneD, Alberta, UG, ALU, FOAM, SP, MMESH, Sub. A grid can be constructed either from parameters or from files. All GridManager classes inherit from the [base GridManager](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/io/grid/gridmanager_base.hh) class.

### Key functionalities

* init():
    - Construct the grid. This needs to be specialized for each grid type.
* grid():
    - Return a reference to the grid.
* loadBalance():
    - Call `loadBalance()` function of the grid.
* makeGridFromFile():
    - Makes a grid from a file. Currently, the formats dgf (Dune grid format), msh (Gmsh mesh format) and vtp/vtu (VTK file format) are supported.
* makeGridFromDgfFile():
    - Makes a grid from a DGF file. This is used by grid managers that only support DGF.
* makeStructuredGrid():
    - Makes a structured cube grid using the structured grid factory of Dune.

### Overview

@mermaid{gridmanager}
