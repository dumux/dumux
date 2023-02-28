## gridview

A gridView is a representation of the grid given a hierarchical view and provides read-only access to certain parts of the grid. A LeafGridView is a view on all elements of the grid without descendants in the hierarchy (which would be the leaves of a hierarchy) while a LevelGridView is a view on all elements of a given level of a refinement hierarchy. The grid view is given by the grid. The leaf grid view can be accessed by `gridManager.grid().leafGridView()`, the level grid view by `gridManager.grid().levelGridView()`. The [default grid view](https://gitlab.dune-project.org/core/dune-grid/-/blob/master/dune/grid/common/defaultgridview.hh) is part of the dune-grid module

### Functionalities

1. grid()
2. size()
3. comm()
4. overlapSize()
5. ghostSize()

### Overview
```{mermaid}
flowchart LR
    A(Grid) -->|"template parameter for declaration"| B(gridManager)
    B -->|"init()"| C(grid)
    C -->|"gridView()"| D(gridView)
    click A "./grid.html"
    click B "./gridmanager.html"
    click C "./grid.html"
    click D "./gridview.html"
```
