# GridView

A gridView is a representation of the grid given a hierarchical view and provides read-only access to certain parts of the grid. A LeafGridView is a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a hierarchy) while a LevelGridView is a view on all elements of a given level of a refinement hierarchy. The grid view is given by the grid. The leaf grid view can be accessed by `gridManager.grid().leafGridView()`, the level grid view by `gridManager.grid().levelGridView()`. The [default grid view](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/defaultgridview.hh) is part of the dune-grid module.

### Key functionalities

* grid():
    - Return a reference to the grid object.
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
