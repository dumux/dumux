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
 * \ingroup BoxDFMModel
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the box discrete fracture model.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_GEOMETRY_HELPER_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_GEOMETRY_HELPER_HH

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

template <class ct>
struct BoxDfmMLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the number of corners in advance (2^(mydim) corners (1<<(mydim))
    // However, on fracture scvs the number might be smaller (use ReservedVector)
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<(mydim)) >;
    };

    // we know all scvfs will have the same geometry type
    template< int mydim >
    struct hasSingleGeometryType
    {
        static const bool v = true;
        static const unsigned int topologyId = Dune::GeometryTypes::cube(mydim).id();
    };
};

//! Create sub control volumes and sub control volume face geometries
template<class GridView, int dim, class ScvType, class ScvfType>
class BoxDfmGeometryHelper;

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class BoxDfmGeometryHelper<GridView, 2, ScvType, ScvfType> : public BoxGeometryHelper<GridView, 2, ScvType, ScvfType>
{
    using ParentType = BoxGeometryHelper<GridView, 2, ScvType, ScvfType>;

    using Intersection = typename GridView::Intersection;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    static constexpr auto dim = GridView::dimension;
    using Scalar = typename GridView::ctype;

public:

    //! Pull up constructor of base class
    using ParentType::ParentType;

    //! Get the corners of the (d-1)-dimensional fracture scvf
    //! The second argument is for compatibility reasons with the 3d case!
    ScvfCornerStorage getFractureScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& isGeom,
                                             unsigned int idxOnIntersection = 0) const
    { return ScvfCornerStorage({isGeom.center()}); }


    //! get fracture scvf normal vector (simply the unit vector of the edge)
    //! The third argument is for compatibility reasons with the 3d case!
    typename ScvfType::Traits::GlobalPosition
    fractureNormal(const ScvfCornerStorage& p,
                   const Intersection& is,
                   unsigned int edgeIndexInIntersection) const
    {
        const auto& geo = this->elementGeometry();
        const auto ref = referenceElement(geo);
        const auto v0 = ref.subEntity(is.indexInInside(), 1, 0, dim);
        const auto v1 = ref.subEntity(is.indexInInside(), 1, 1, dim);
        auto normal = geo.corner(v1) - geo.corner(v0);
        normal /= normal.two_norm();
        return normal;
    }
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class BoxDfmGeometryHelper<GridView, 3, ScvType, ScvfType> : public BoxGeometryHelper<GridView, 3, ScvType, ScvfType>
{
    using ParentType = BoxGeometryHelper<GridView, 3, ScvType, ScvfType>;

    using Intersection = typename GridView::Intersection;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;
    using Scalar = typename GridView::ctype;

public:

    //! Pull up constructor of base class
    using ParentType::ParentType;

    //! Create the sub control volume face geometries on an intersection marked as fracture
    ScvfCornerStorage getFractureScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& isGeom,
                                             unsigned int edgeIndexInIntersection) const
    {
        const auto& geo = this->elementGeometry();
        const auto refElement = referenceElement(geo);
        const auto faceRefElem = referenceElement(isGeom);

        // create point vector for this geometry
        typename ScvfType::Traits::GlobalPosition pi[9];

        // the facet center
        pi[0] = isGeom.center();

        // facet edge midpoints
        const auto idxInInside = is.indexInInside();
        for (int i = 0; i < faceRefElem.size(1); ++i)
        {
            const auto edgeIdxLocal = refElement.subEntity(idxInInside, 1, i, dim-1);
            pi[i+1] = geo.global(refElement.position(edgeIdxLocal, dim-1));
        }

        // proceed according to number of corners
        const auto corners = isGeom.corners();
        switch (corners)
        {
            case 3: // triangle
            {
                //! Only build the maps the first time we encounter a triangle
                static const std::uint8_t fo = 1; //!< face offset in point vector p
                static const std::uint8_t map[3][2] =
                {
                    {fo+0, 0},
                    {0, fo+1},
                    {fo+2, 0}
                };

                return ScvfCornerStorage{ {pi[map[edgeIndexInIntersection][0]],
                                           pi[map[edgeIndexInIntersection][1]]} };
            }
            case 4: // quadrilateral
            {
                //! Only build the maps the first time we encounter a quadrilateral
                static const std::uint8_t fo = 1; //!< face offset in point vector p
                static const std::uint8_t map[4][2] =
                {
                    {0, fo+0},
                    {fo+1, 0},
                    {fo+2, 0},
                    {0, fo+3}
                };

                return ScvfCornerStorage{ {pi[map[edgeIndexInIntersection][0]],
                                           pi[map[edgeIndexInIntersection][1]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box fracture scvf geometries for dim=" << dim
                                                                << " dimWorld=" << dimWorld
                                                                << " corners=" << corners);
        }
    }

    //! get fracture scvf normal vector
    typename ScvfType::Traits::GlobalPosition
    fractureNormal(const ScvfCornerStorage& scvfCorners,
                   const Intersection& is,
                   unsigned int edgeIndexInIntersection) const
    {
        const auto& geo = this->elementGeometry();
        const auto refElement = referenceElement(geo);

        // first get the intersection corners (maximum "4" is for quadrilateral face)
        typename ScvfType::Traits::GlobalPosition c[4];

        const auto corners = refElement.size(is.indexInInside(), 1, dim);
        for (int i = 0; i < corners; ++i)
        {
            const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, i, dim);
            c[i] = geo.corner(vIdxLocal);
        }

        // compute edge vector depending on number of corners
        const auto gridEdge = [&] ()
        {
            // triangles
            if (corners == 3)
            {
                if (edgeIndexInIntersection == 0) return c[1]-c[0];
                else if (edgeIndexInIntersection == 1) return c[2]-c[0];
                else if (edgeIndexInIntersection == 2) return c[2]-c[1];
                else DUNE_THROW(Dune::InvalidStateException, "Invalid edge index");
            }
            else if (corners == 4)
            {
                if (edgeIndexInIntersection == 0) return c[2]-c[0];
                else if (edgeIndexInIntersection == 1) return c[3]-c[1];
                else if (edgeIndexInIntersection == 2) return c[1]-c[0];
                else if (edgeIndexInIntersection == 3) return c[3]-c[2];
                else DUNE_THROW(Dune::InvalidStateException, "Invalid edge index");
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid face geometry");
        } ();

        // compute lower edge of the scvf
        assert(scvfCorners.size() == 2);
        const auto scvfEdge = scvfCorners[1]-scvfCorners[0];

        // compute scvf normal via 2 cross products
        const auto faceN = crossProduct(gridEdge, scvfEdge);
        auto n = crossProduct(scvfEdge, faceN);
        n /= n.two_norm();
        return n;
    }
};

} // end namespace Dumux

#endif
