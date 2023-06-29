# GridView

A *grid view* is a representation of the grid that allows for read-only access to certain parts of the *grid* from which it is obtained. A *leaf grid view* is a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a grid hierarchy) while a *level grid view* is a view on all elements of a given level of a refinement hierarchy. Refinements can be either adaptive (for grid adaptvive methods) or static (e.g. for multigrid methods). Grid view objects are accessed via the grid object as follows: The leaf grid view can be accessed by `gridManager.grid().leafGridView()`, the level grid view by `gridManager.grid().levelGridView()`. The [default grid view](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/defaultgridview.hh) is part of the dune-grid module. Note that grid views provide *forward iteration* acess only, so random access operators are not provided.

### Key functionalities

* grid():
    - Return a reference to the grid object from which the `GridView` object is obtained.
* size(...):
    - Return number of grid or leaf entities of a given `codim` on a given `level`. Also see the explanation of the function `size()` in @ref grid.
* comm():
    - Return the communication object which manages parallelization in case of a parallel grid.
* overlapSize(codim):
    - Return the size of the overlap region for a given `codim` and grid view.
* ghostSize(codim):
    - Return the size of the ghost region for a given `codim` and grid view.

### Overview

@mermaid{gridview}
