# GridGeometry

A grid geometry is an abstraction of a finite-volume discretization and differs from discretization to discretization. Given a @ref gridview "grid view", the grid geometry constructs all geometrical and topological data necessary to evaluate the discrete equation resulting from a given finite volume scheme. All grid geometries inherit from the [basic grid geometry](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/basicgridgeometry.hh).

### Key functionalities

* update(gridView):
    - Update internal state after grid changed.
* gridView():
    - Return the gridView this grid geometry object lives on.
* vertexMapper():
    - Returns the mapper for vertices to indices for possibly adaptive grids.
* elementMapper():
    - Returns the mapper for elements to indices for possibly adaptive grids.
* boundingBoxTree():
    - Returns the bounding box tree of the grid. This is a hierarchy of axis-aligned boxes, on the lowest level of the hierarchie the bounding boxes of the individual elements can be found while on the highest level (root) one single box bounding all elements can be found.. For example, this tree can be used to conduct binary searches in order to find the element which is located at a specific position.
* elementMap():
    - Returns the element index to element map.
* element(elementIndex):
    - Get an element from a global `elementIndex`.
* bBoxMin():
    - The coordinate of the corner of the `GridView's` bounding box with the smallest values. For a cartesian grid, this would be the coordinate of the lower-left corner.
* bBoxMax():
    - The coordinate of the corner of the `GridView's` bounding box with the largest values. For a cartesian grid, this would be the coordinate of the upper-right corner.

### Overview

@mermaid{gridgeometry}