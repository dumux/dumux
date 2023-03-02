## GridGeometry {#gridgeometry}
<!-- @page gridGeometry GridGeometry -->


A grid geometry is an abstraction of a finite-volume discretization. Given a @ref gridView, the grid geometry constructs all geometrical and topological data necessary to evaluate the discrete equation resulting from a given finite volume scheme.
Since Dumux uses elementwise assembly one needs to be able to access the geometry from a element perspective.
This is implemented using the @ref fvElementGeometry` object.

### Functionalities

1. elementMapper()
2. vertexMapper()
3. boundingBox()
