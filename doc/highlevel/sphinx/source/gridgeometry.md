## gridgeometry

A grid geometry is an abstraction of a finite-volume discretization and differs from discretization to discretization. Given a {ref}`gridview`, the grid geometry constructs all geometrical and topological data necessary to evaluate the discrete equation resulting from a given finite volume scheme. All grid geometries inherit from the [basic grid geometry](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/basicgridgeometry.hh).

### Functionalities

1. update()
2. gridView()
3. vertexMapper()
4. elementMapper()
5. boundingBoxTree()
6. elementMap()
7. element()
8. bBoxMin()
9. bBoxMax()

### Overview
```{mermaid}
flowchart LR
    A(Grid) -->|"template parameter for declaration"| B(gridManager)
    B -->|"init()"| C(grid)
    C -->|"gridView()"| D(gridView)
    D --> |"argument for constructor"| E(gridgeometry)
    click A "./grid.html"
    click B "./gridmanager.html"
    click C "./grid.html"
    click D "./gridview.html"
    click E "./gridgeometry.html"
```
