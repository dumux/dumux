# GridGeometry

A grid geometry is an abstraction of a finite-volume discretization and differs fora each discretization method. Given a @ref gridview "grid view", the grid geometry constructs all geometrical and topological data necessary to evaluate the discrete equations resulting from a given finite volume scheme. All grid geometries inherit from the [basic grid geometry](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/basicgridgeometry.hh).

### Example

```cpp
...
GridManager<Grid> gridManager;
gridManager.init();
auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView()); // We, create the finite volume grid geometry from the (leaf) grid view.
...
```

In this case, a shared pointer to the `gridGeometry` object is created, given the `leafGridView` of the `grid`. Also, the `gridGeometry` is used as an argument for initialising objects from other classes, such as a `problem` object or a `gridVariables` object.

### Key functionalities

* See Dumux::BasicGridGeometry.

### Overview

@mermaid{gridgeometry}
