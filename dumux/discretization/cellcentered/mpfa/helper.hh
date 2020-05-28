// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to get data required for mpfa scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Dimension-specific helper class to get data required for mpfa scheme.
 */
template<class GridGeometry, int dim, int dimWorld>
class MpfaDimensionHelper;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief  Dimension-specific mpfa helper class for dim == 2 & dimWorld == 2
 */
template<class GridGeometry>
class MpfaDimensionHelper<GridGeometry, /*dim*/2, /*dimWorld*/2>
{
    using GridView = typename GridGeometry::GridView;
    using CoordScalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfCornerVector = typename SubControlVolumeFace::Traits::CornerStorage;

    // Container to store the positions of intersections required for scvf
    // corner computation. In 2d, these are the center plus the two corners
    using ScvfPositionsOnIntersection = std::array<GlobalPosition, 3>;

public:
    /*!
     * \brief Calculates the inner normal vectors to a given scv basis.
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static ScvBasis calculateInnerNormals(const ScvBasis& scvBasis)
    {
        static const Dune::FieldMatrix<CoordScalar, 2, 2> R = {{0.0, 1.0}, {-1.0, 0.0}};

        ScvBasis innerNormals;
        R.mv(scvBasis[1], innerNormals[0]);
        R.mv(scvBasis[0], innerNormals[1]);

        // adjust sign depending on basis being a RHS
        if (!isRightHandSystem(scvBasis))
            innerNormals[0] *= -1.0;
        else
            innerNormals[1] *= -1.0;

        return innerNormals;
    }

    /*!
     * \brief Calculates the determinant of an scv basis.
     *        This is equal to the cross product for dim = dimWorld = 2
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static CoordScalar calculateDetX(const ScvBasis& scvBasis)
    {
        using std::abs;
        return abs(crossProduct<CoordScalar>(scvBasis[0], scvBasis[1]));
    }

    /*!
     * \brief Returns the global number of scvfs in the grid. Assumes the grid to be made up of only
     *        basic geometry types. Overlad this function if you want to use different geometry types.
     * \param gridView The grid view to be checked
     */
    static std::size_t getGlobalNumScvf(const GridView& gridView)
    {
        assert(gridView.size(Dune::GeometryTypes::triangle)
               + gridView.size(Dune::GeometryTypes::quadrilateral) == gridView.size(0));

        return gridView.size(Dune::GeometryTypes::triangle)
                 * getNumLocalScvfs(Dune::GeometryTypes::triangle)
               + gridView.size(Dune::GeometryTypes::quadrilateral)
                 * getNumLocalScvfs(Dune::GeometryTypes::quadrilateral);
    }

    /*!
     * \brief Checks whether or not a given scv basis forms a right hand system.
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static bool isRightHandSystem(const ScvBasis& scvBasis)
    { return !std::signbit(crossProduct<CoordScalar>(scvBasis[0], scvBasis[1])); }

    /*!
     * \brief Returns a vector containing the positions on a given intersection
     *        that are relevant for scvf corner computation.
     *        Ordering -> 1: facet center, 2: the two facet corners
     *
     * \param eg Geometry of the element the facet is embedded in
     * \param refElement Reference element of the element the facet is embedded in
     * \param indexInElement The local index of the facet in the element
     * \param numCorners The number of corners on the facet
     */
    template<class ElementGeometry, class ReferenceElement>
    static ScvfPositionsOnIntersection computeScvfCornersOnIntersection(const ElementGeometry& eg,
                                                                        const ReferenceElement& refElement,
                                                                        unsigned int indexInElement,
                                                                        unsigned int numCorners)
    {
        ScvfPositionsOnIntersection p;

        // compute facet center and corners
        p[0] = 0.0;
        for (unsigned int c = 0; c < numCorners; ++c)
        {
            // codim = 1, dim = 2
            p[c+1] = eg.global(refElement.position(refElement.subEntity(indexInElement, 1, c, 2), 2));
            p[0] += p[c+1];
        }
        p[0] /= numCorners;

        return p;
    }

    /*!
     * \brief Returns the corners of the sub control volume face constructed
     *        in a corner (vertex) of an intersection.
     *
     * \param p Container with all scvf corners of the intersection
     * \param numIntersectionCorners Number of corners of the intersection (required in 3d)
     * \param cornerIdx Local vertex index on the intersection
     */
    static ScvfCornerVector getScvfCorners(const ScvfPositionsOnIntersection& p,
                                           unsigned int numIntersectionCorners,
                                           unsigned int cornerIdx)
    {
        // make sure the given input is admissible
        assert(cornerIdx < 2 && "provided index exceeds the number of corners of facets in 2d");

        // create & return the scvf corner vector
        return cornerIdx == 0 ? ScvfCornerVector({p[0], p[1]})
                              : ScvfCornerVector({p[0], p[2]});
    }

    /*!
     * \brief Calculates the area of an scvf.
     * \param scvfCorners Container with the corners of the scvf
     */
    static CoordScalar computeScvfArea(const ScvfCornerVector& scvfCorners)
    { return (scvfCorners[1]-scvfCorners[0]).two_norm(); }

    /*!
     * \brief Calculates the number of scvfs in a given element geometry type.
     * \param gt The element geometry type
     */
    static std::size_t getNumLocalScvfs(const Dune::GeometryType& gt)
    {
        if (gt == Dune::GeometryTypes::triangle)
            return 6;
        else if (gt == Dune::GeometryTypes::quadrilateral)
            return 8;
        else
            DUNE_THROW(Dune::NotImplemented, "Mpfa for 2d geometry type " << gt);
    }
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief  Dimension-specific mpfa helper class for dim == 2 & dimWorld == 3.
 *         Reuses some functionality of the specialization for dim = dimWorld = 2
 */
template<class GridGeometry>
class MpfaDimensionHelper<GridGeometry, /*dim*/2, /*dimWorld*/3>
: public MpfaDimensionHelper<GridGeometry, 2, 2>
{
    using GridView = typename GridGeometry::GridView;
    using CoordScalar = typename GridView::ctype;
public:

    /*!
     * \brief Calculates the inner normal vectors to a given scv basis.
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static ScvBasis calculateInnerNormals(const ScvBasis& scvBasis)
    {
        // compute vector normal to the basis plane
        const auto normal = [&scvBasis] ()
        {
            auto n = crossProduct<CoordScalar>(scvBasis[0], scvBasis[1]);
            n /= n.two_norm();
            return n;
        } ();

        // compute inner normals using the normal vector
        ScvBasis innerNormals;
        innerNormals[0] = crossProduct<CoordScalar>(scvBasis[1], normal);
        innerNormals[1] = crossProduct<CoordScalar>(normal, scvBasis[0]);

        return innerNormals;
    }

    /*!
     * \brief Calculates the determinant of an scv basis.
     *        For dim = 2 < dimWorld = 3 this is actually not the determinant of the
     *        basis but it is simply the area of the parallelofram spanned by the
     *        basis vectors.
     *
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static CoordScalar calculateDetX(const ScvBasis& scvBasis)
    {
        using std::abs;
        return abs(crossProduct<CoordScalar>(scvBasis[0], scvBasis[1]).two_norm());
    }

    /*!
     * \brief Checks whether or not a given scv basis forms a right hand system.
     *        Note that for dim = 2 < dimWorld = 3 the bases forming a right hand system
     *        are not unique.
     *
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static constexpr bool isRightHandSystem(const ScvBasis& scvBasis) { return true; }
};
/*!
 * \ingroup CCMpfaDiscretization
 * \brief   Dimension-specific mpfa helper class for dim == 3 & dimWorld == 3.
 *
 */
template<class GridGeometry>
class MpfaDimensionHelper<GridGeometry, /*dim*/3, /*dimWorld*/3>
{
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfCornerVector = typename SubControlVolumeFace::Traits::CornerStorage;

    // Be picky about the dimensions
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // container to store the positions of intersections required for
    // scvf corner computation. Maximum number of points needed is 9
    // for the supported geometry types (quadrilateral facet)
    using ScvfPositionsOnIntersection = std::array<GlobalPosition, 9>;

public:

    /*!
     * \brief Calculates the inner normal vectors to a given scv basis.
     *
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static ScvBasis calculateInnerNormals(const ScvBasis& scvBasis)
    {
        ScvBasis innerNormals;

        innerNormals[0] = crossProduct<CoordScalar>(scvBasis[1], scvBasis[2]);
        innerNormals[1] = crossProduct<CoordScalar>(scvBasis[2], scvBasis[0]);
        innerNormals[2] = crossProduct<CoordScalar>(scvBasis[0], scvBasis[1]);

        if (!isRightHandSystem(scvBasis))
            std::for_each(innerNormals.begin(), innerNormals.end(), [] (auto& n) { n *= -1.0; });

        return innerNormals;
    }

    /*!
     * \brief Calculates the determinant of an scv basis.
     *        This is equal to the cross product for dim = dimWorld = 2
     *
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static CoordScalar calculateDetX(const ScvBasis& scvBasis)
    {
        using std::abs;
        return abs(tripleProduct<CoordScalar>(scvBasis[0], scvBasis[1], scvBasis[2]));
    }

    /*!
     * \brief Returns the global number of scvfs in the grid. Assumes the grid to be made up of only
     *        basic geometry types. Overload this function if you want to use different geometry types.
     *
     * \param gridView The grid view to be checked
     */
    static std::size_t getGlobalNumScvf(const GridView& gridView)
    {
        assert(gridView.size(Dune::GeometryTypes::tetrahedron)
               + gridView.size(Dune::GeometryTypes::pyramid)
               + gridView.size(Dune::GeometryTypes::prism)
               + gridView.size(Dune::GeometryTypes::hexahedron) == gridView.size(0));

        return gridView.size(Dune::GeometryTypes::tetrahedron)
                 * getNumLocalScvfs(Dune::GeometryTypes::tetrahedron)
               + gridView.size(Dune::GeometryTypes::pyramid)
                 * getNumLocalScvfs(Dune::GeometryTypes::pyramid)
               + gridView.size(Dune::GeometryTypes::prism)
                 * getNumLocalScvfs(Dune::GeometryTypes::prism)
               + gridView.size(Dune::GeometryTypes::hexahedron)
                 * getNumLocalScvfs(Dune::GeometryTypes::hexahedron);
    }

    /*!
     * \brief Checks whether or not a given scv basis forms a right hand system.
     * \param scvBasis The basis of an scv
     */
    template< class ScvBasis >
    static bool isRightHandSystem(const ScvBasis& scvBasis)
    { return !std::signbit(tripleProduct<CoordScalar>(scvBasis[0], scvBasis[1], scvBasis[2])); }

    /*!
     * \brief Returns a vector containing the positions on a given intersection
     *        that are relevant for scvf corner computation.
     *        Ordering -> 1: facet center, 2: facet corners, 3: edge centers
     *
     * \param eg Geometry of the element the facet is embedded in
     * \param refElement Reference element of the element the facet is embedded in
     * \param indexInElement The local index of the facet in the element
     * \param numCorners The number of corners on the facet
     */
    template<class ElementGeometry, class ReferenceElement>
    static ScvfPositionsOnIntersection computeScvfCornersOnIntersection(const ElementGeometry& eg,
                                                                        const ReferenceElement& refElement,
                                                                        unsigned int indexInElement,
                                                                        unsigned int numCorners)
    {
        // The size of ScvfPositionsOnIntersection doesn't allow for faces with more than four corners!
        ScvfPositionsOnIntersection p;
        if (numCorners > 4)
            DUNE_THROW(Dune::InvalidStateException, "Mpfa implementation cannot handle faces with more than 4 corners");

        // compute facet center and corners
        p[0] = 0.0;
        for (unsigned int c = 0; c < numCorners; ++c)
        {
            // codim = 1, dim = 3
            p[c+1] = eg.global(refElement.position(refElement.subEntity(indexInElement, 1, c, 3), 3));
            p[0] += p[c+1];
        }
        p[0] /= numCorners;

        // proceed according to number of corners
        switch (numCorners)
        {
            case 3: // triangle
            {
                // add edge midpoints, there are 3 for triangles
                p[numCorners+1] = p[2] + p[1];
                p[numCorners+1] /= 2;
                p[numCorners+2] = p[3] + p[1];
                p[numCorners+2] /= 2;
                p[numCorners+3] = p[3] + p[2];
                p[numCorners+3] /= 2;
                return p;
            }
            case 4: // quadrilateral
            {
                // add edge midpoints, there are 4 for quadrilaterals
                p[numCorners+1] = p[3] + p[1];
                p[numCorners+1] /= 2;
                p[numCorners+2] = p[4] + p[2];
                p[numCorners+2] /= 2;
                p[numCorners+3] = p[2] + p[1];
                p[numCorners+3] /= 2;
                p[numCorners+4] = p[4] + p[3];
                p[numCorners+4] /= 2;
                return p;
            }
            default:
                DUNE_THROW(Dune::NotImplemented,
                           "Mpfa scvf corners for dim = 3, dimWorld = 3, corners = " << numCorners);
        }
    }

    /*!
     * \brief Returns the corners of the sub control volume face constructed
     *        in a corner (vertex) of an intersection.
     *
     * \param p Container with all scvf corners of the intersection
     * \param numIntersectionCorners Number of corners of the intersection
     * \param cornerIdx Local vertex index on the intersection
     */
    static ScvfCornerVector getScvfCorners(const ScvfPositionsOnIntersection& p,
                                           unsigned int numIntersectionCorners,
                                           unsigned int cornerIdx)
    {
        // proceed according to number of corners
        // we assume the ordering according to the above method computeScvfCornersOnIntersection()
        switch (numIntersectionCorners)
        {
            case 3: // triangle
            {
                // Only build the maps the first time we encounter a triangle
                static const std::uint8_t vo = 1; // vertex offset in point vector p
                static const std::uint8_t eo = 4; // edge offset in point vector p
                static const std::uint8_t map[3][4] =
                {
                    {0, eo+1, eo+0, vo+0},
                    {0, eo+0, eo+2, vo+1},
                    {0, eo+2, eo+1, vo+2}
                };

                return ScvfCornerVector( {p[map[cornerIdx][0]],
                                          p[map[cornerIdx][1]],
                                          p[map[cornerIdx][2]],
                                          p[map[cornerIdx][3]]} );
            }
            case 4: // quadrilateral
            {
                // Only build the maps the first time we encounter a quadrilateral
                static const std::uint8_t vo = 1; // vertex offset in point vector p
                static const std::uint8_t eo = 5; // face offset in point vector p
                static const std::uint8_t map[4][4] =
                {
                    {0, eo+0, eo+2, vo+0},
                    {0, eo+2, eo+1, vo+1},
                    {0, eo+3, eo+0, vo+2},
                    {0, eo+1, eo+3, vo+3}
                };

                return ScvfCornerVector( {p[map[cornerIdx][0]],
                                          p[map[cornerIdx][1]],
                                          p[map[cornerIdx][2]],
                                          p[map[cornerIdx][3]]} );
            }
            default:
                DUNE_THROW(Dune::NotImplemented,
                           "Mpfa scvf corners for dim = 3, dimWorld = 3, corners = " << numIntersectionCorners);
        }
    }

    /*!
     * \brief Calculates the area of an scvf.
     * \param scvfCorners Container with the corners of the scvf
     */
    static CoordScalar computeScvfArea(const ScvfCornerVector& scvfCorners)
    {
        // after Wolfram alpha quadrilateral area
        return 0.5*Dumux::crossProduct(scvfCorners[3]-scvfCorners[0], scvfCorners[2]-scvfCorners[1]).two_norm();
    }

    /*!
     * \brief Calculates the number of scvfs in a given element geometry type.
     *
     * \param gt The element geometry type
     */
    static std::size_t getNumLocalScvfs(const Dune::GeometryType& gt)
    {
        if (gt == Dune::GeometryTypes::tetrahedron)
            return 12;
        else if (gt == Dune::GeometryTypes::pyramid)
            return 16;
        else if (gt == Dune::GeometryTypes::prism)
            return 18;
        else if (gt == Dune::GeometryTypes::hexahedron)
            return 24;
        else
            DUNE_THROW(Dune::NotImplemented, "Mpfa for 3d geometry type " << gt);
    }
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to get the required information on an interaction volume.
 *
 * \tparam GridGeometry The finite volume grid geometry
 */
template<class GridGeometry>
class CCMpfaHelper : public MpfaDimensionHelper<GridGeometry,
                                                GridGeometry::GridView::dimension,
                                                GridGeometry::GridView::dimensionworld>
{
    using PrimaryInteractionVolume = typename GridGeometry::GridIVIndexSets::PrimaryInteractionVolume;
    using SecondaryInteractionVolume = typename GridGeometry::GridIVIndexSets::SecondaryInteractionVolume;

    using VertexMapper = typename GridGeometry::VertexMapper;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ScvfCornerVector = typename SubControlVolumeFace::Traits::CornerStorage;

    using GridView = typename GridGeometry::GridView;
    static constexpr int dim = GridView::dimension;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CoordScalar = typename GridView::ctype;

public:
    /*!
     * \brief Calculates the integration point on an scvf.
     *
     * \param scvfCorners Container with the corners of the scvf
     * \param q Parameterization of the integration point on the scvf
     */
    template< class Scalar >
    static GlobalPosition getScvfIntegrationPoint(const ScvfCornerVector& scvfCorners, Scalar q)
    {
        // ordering -> first corner: facet center, last corner: vertex
        if (q == 0.0)
            return scvfCorners[0];

        auto ip = scvfCorners.back() - scvfCorners.front();
        ip *= q;
        ip += scvfCorners[0];
        return ip;
    }

    /*!
     * \brief Returns a vector which maps true to each vertex on processor boundaries and false otherwise
     * \todo TODO: The name of this function is not so good, as these are not ghost vertices according
     *             to the Dune definition of ghost entities. Moreover, it should be tried to make MPFA work
     *             also with ghost entities.
     */
    static std::vector<bool> findGhostVertices(const GridView& gridView, const VertexMapper& vertexMapper)
    {
        std::vector<bool> ghostVertices(gridView.size(dim), false);

        // if not run in parallel, skip the rest
        if (Dune::MPIHelper::getCollectiveCommunication().size() == 1)
            return ghostVertices;

        // mpfa methods cannot yet handle ghost cells
        if (gridView.ghostSize(0) > 0)
            DUNE_THROW(Dune::InvalidStateException, "Mpfa methods in parallel do not work with ghost cells. Use overlap cells instead.");

        // mpfa methods have to have overlapping cells
        if (gridView.overlapSize(0) == 0)
            DUNE_THROW(Dune::InvalidStateException, "Grid no overlaping cells. This is required by mpfa methods in parallel.");

        for (const auto& element : elements(gridView))
        {
            for (const auto& is : intersections(gridView, element))
            {
                if (!is.neighbor() && !is.boundary())
                {
                    const auto refElement = referenceElement(element);

                    for (int isVertex = 0; isVertex < is.geometry().corners(); ++isVertex)
                    {
                        const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, isVertex, dim);
                        const auto vIdxGlobal = vertexMapper.subIndex(element, vIdxLocal, dim);
                        ghostVertices[vIdxGlobal] = true;
                    }
                }
            }
        }

        return ghostVertices;
    }

    //! Returns whether or not secondary interaction volumes have to be considered in the model.
    //! This is always the case when the specified types for the interaction volumes differ.
    static constexpr bool considerSecondaryIVs()
    { return !std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value; }

    //! returns whether or not a value exists in a vector
    template<typename V, typename T>
    static bool vectorContainsValue(const V& vector, const T value)
    { return std::find(vector.begin(), vector.end(), value) != vector.end(); }
};

} // end namespace Dumux

#endif
