@addtogroup Geometry

The geometry algorithms find the entities of a grid, or of any other set of geometries, that intersect
a point, a geometry or the entities of a second set, and they compute these intersections.
Point sources use them to find the elements containing a source,
and the embedded and mixed-dimensional couplings (@ref MultiDomain) use them to compute the coupling interfaces.

# Overview {#geometry-overview}

The intersection algorithms are built up in layers:

| Layer | Header | Main classes and functions |
|:------|:-------|:---------------------------|
| entity sets | `dumux/geometry/geometricentityset.hh` | @ref Dumux::GridViewGeometricEntitySet, @ref Dumux::GeometriesEntitySet, @ref Dumux::FixedSizeGeometriesEntitySet, @ref Dumux::SingleGeometryEntitySet |
| bounding box tree | `dumux/geometry/boundingboxtree.hh` | @ref Dumux::BoundingBoxTree |
| queries on trees | `dumux/geometry/intersectingentities.hh` | `Dumux::intersectingEntities`, `Dumux::intersectingBoxes`, @ref Dumux::IntersectionInfo |
| intersection of two geometries | `dumux/geometry/geometryintersection.hh` | @ref Dumux::GeometryIntersection, the policies in `Dumux::IntersectionPolicy` |
| triangulation of intersections | `dumux/geometry/triangulation.hh` | `Dumux::triangulate` |
| intersections of two entity sets | `dumux/geometry/intersectionentityset.hh` | @ref Dumux::IntersectionEntitySet |
| intersections of two grids | `dumux/multidomain/glue.hh` | `Dumux::makeGlue`, used by the projections in `dumux/discretization/projection/projector.hh` |

The distributed (MPI-parallel) counterparts are @ref Dumux::DistributedBoundingBoxTree,
the `intersectingEntities` overloads in `dumux/geometry/distributedintersectingentities.hh`
and @ref Dumux::DistributedIntersectionEntitySet.

# Entity sets and bounding box trees {#geometry-entity-sets-and-trees}

An entity set gives access to a set of geometries through consecutive indices.
@ref Dumux::GridViewGeometricEntitySet wraps the entities of a grid view
(by default its elements, with the indices of a Dune mapper),
@ref Dumux::GeometriesEntitySet stores a vector of geometries.
The queries below act on a @ref Dumux::BoundingBoxTree built for an entity set.

@snippet test_intersectingentities.cc grid tree

@snippet test_intersectingentities.cc geometries tree

If the geometries of the set change, the tree has to be built again with `build(entitySet)`.
The grid geometries of the discretizations provide a tree of the elements of their grid view,
`gridGeometry.boundingBoxTree()`.

# Queries {#geometry-queries}

The function `intersectingEntities` finds the entities of a tree that intersect a point,
a geometry or the entities of a second tree:

| Query | Result |
|:------|:-------|
| `intersectingEntities(point, tree)` | the indices of all entities containing the point |
| `intersectingEntities(geometry, tree)` | an @ref Dumux::IntersectionInfo per intersection, with the entity index in `second()` |
| `intersectingEntities(treeA, treeB)` | an @ref Dumux::IntersectionInfo per intersection, with the entity indices in `first()` and `second()` |
| `intersectingBoxes(treeA, treeB)` | the index pairs of all entities whose bounding boxes overlap (no geometry test) |
| `intersectingEntityCartesianGrid(point, min, max, cells)` | the index of the cell of a structured grid that contains a point inside the grid (no tree) |

The point query returns all entities that contain the point, so a point on the boundary between
entities is found in each of them:

@snippet test_intersectingentities.cc point query

An @ref Dumux::IntersectionInfo holds the corners of the intersection geometry, `corners()`.
Intersections of polygons and polyhedra are split into triangles or tetrahedra,
and every simplex is returned as one @ref Dumux::IntersectionInfo,
so the same pair of entities can occur several times:

@snippet test_intersectingentities.cc tree query

The dimension of the computed intersections is given by an intersection policy (see below),
passed as the last argument. By default, it is the smaller dimension of the two geometries,
e.g. segments for a polyline and the cells of a grid:

@snippet test_intersectingentities.cc policy

A policy other than the default is useful, e.g., for two sets of segments: the default policy computes
the segments where they overlap, `IntersectionPolicy::PointPolicy` the points where they cross.

# Intersection of two geometries {#geometry-intersection-of-two-geometries}

@ref Dumux::GeometryIntersection tests and intersects two geometries with the Dune geometry interface.
The intersection policy is the third template argument:

| Policy | `Intersection` type | Content |
|:-------|:--------------------|:--------|
| `IntersectionPolicy::PointPolicy` | `Dune::FieldVector` | the intersection point |
| `IntersectionPolicy::SegmentPolicy` | `std::array` of two points | the end points of the intersection segment |
| `IntersectionPolicy::PolygonPolicy` | `std::vector` of points | the corners of the convex intersection polygon, ordered along its boundary |
| `IntersectionPolicy::PolyhedronPolicy` | `std::vector` of points | the corners of the convex intersection polyhedron, unordered (the intersection is their convex hull) |

@snippet test_intersectingentities.cc geometry intersection

The following combinations are implemented, for both orders of the two geometries.
Other combinations throw `Dune::NotImplemented`.

| World dimension | Geometries | Point | Segment | Polygon | Polyhedron |
|:---:|:-----------------------|:---:|:---:|:---:|:---:|
| 2 | segment, segment       | yes | yes |     |     |
| 2 | polygon, segment       | yes | yes |     |     |
| 2 | polygon, polygon       |     |     | yes |     |
| 3 | segment, segment       | yes | yes |     |     |
| 3 | polygon, segment       | yes | yes |     |     |
| 3 | polygon, polygon       |     |     | yes (coplanar polygons) |     |
| 3 | polyhedron, segment    | yes | yes |     |     |
| 3 | polyhedron, polygon    |     |     | yes |     |
| 3 | polyhedron, polyhedron |     |     |     | yes |

A point-like intersection is a single point: segments that overlap a segment or a polygon edge along a part
of their length have no point-like intersection with it.

The algorithms assume convex geometries with straight edges and planar faces:
segments, triangles and quadrilaterals, tetrahedra, pyramids, prisms and hexahedra.
Their corners may be numbered in either orientation (e.g. mirrored elements).
Two neighboring elements that only share a face have no polygon- or polyhedron-like intersection.

Without a tree, `intersectsPointGeometry` (`dumux/geometry/intersectspointgeometry.hh`)
tests whether a point lies in a geometry.

# Tolerances {#geometry-tolerances}

All tests use tolerances relative to the size of the involved geometries.
A point outside a geometry by less than this tolerance counts as inside,
and points closer to each other than this tolerance count as one point.

# Intersections of two entity sets and of two grids {#geometry-entity-set-intersections}

@ref Dumux::IntersectionEntitySet holds the intersections of two entity sets as entities,
which have a geometry and know the intersecting entities of both sets:

```cpp
IntersectionEntitySet<DomainEntitySet, TargetEntitySet> glue;
glue.build(domainTree, targetTree);
for (const auto& is : intersections(glue))
{
    const auto geometry = is.geometry();
    const auto domainEntity = is.domainEntity();
    const auto targetEntity = is.targetEntity();
    // ...
}
```

If the two sets have different dimensions, a geometrically identical intersection of a lower-dimensional
entity with several higher-dimensional entities (e.g. a fracture segment on the face between two elements)
is stored once, with all of them as neighbors.
`Dumux::makeGlue(domainGridGeometry, targetGridGeometry)` (`dumux/multidomain/glue.hh`)
builds such a set for the grids of two grid geometries.

# Further geometric functions {#geometry-further-functions}

| Header | Content |
|:-------|:--------|
| `dumux/geometry/distance.hh` | distances between points, lines, segments, triangles, polygons and geometries |
| `dumux/geometry/distancefield.hh` | @ref Dumux::AABBDistanceField, distance from a point to the closest of a set of geometries |
| `dumux/geometry/volume.hh` | `convexPolytopeVolume`, `volume` of a geometry |
| `dumux/geometry/center.hh`, `diameter.hh`, `normal.hh` | center of a set of corners, diameter of a geometry, normal vectors |
| `dumux/geometry/boundingsphere.hh`, `circumsphere.hh` | bounding sphere of a geometry or a point cloud, circumsphere of a triangle |
| `dumux/geometry/grahamconvexhull.hh` | convex hull of a set of points in a plane |
| `dumux/geometry/makegeometry.hh` | quadrilaterals from four points in arbitrary order |
| `dumux/geometry/refinementquadraturerule.hh` | quadrature rule on a virtually refined geometry |
