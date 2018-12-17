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
 * \ingroup BoxDiscretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the box scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct BoxDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar>>;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDefaultScvfGeometryTraits<GV> >
class BoxSubControlVolumeFace
: public SubControlVolumeFaceBase<BoxSubControlVolumeFace<GV, T>, T>
{
    using ThisType = BoxSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class GeometryHelper, class Element>
    BoxSubControlVolumeFace(const GeometryHelper& geometryHelper,
                            const Element& element,
                            const typename Element::Geometry& elemGeometry,
                            GridIndexType scvfIndex,
                            std::vector<LocalIndexType>&& scvIndices,
                            bool boundary = false)
    : corners_(geometryHelper.getScvfCorners(scvfIndex)),
      center_(0.0),
      unitOuterNormal_(geometryHelper.normal(corners_, scvIndices)),
      area_(geometryHelper.scvfArea(corners_)),
      scvfIndex_(scvfIndex),
      scvIndices_(std::move(scvIndices)),
      boundary_(boundary)
    , boundaryFlag_{}
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for boundary scvfs
    template<class GeometryHelper, class Intersection>
    BoxSubControlVolumeFace(const GeometryHelper& geometryHelper,
                            const Intersection& intersection,
                            const typename Intersection::Geometry& isGeometry,
                            LocalIndexType indexInIntersection,
                            GridIndexType scvfIndex,
                            std::vector<LocalIndexType>&& scvIndices,
                            bool boundary = false)
    : corners_(geometryHelper.getBoundaryScvfCorners(intersection, isGeometry, indexInIntersection)),
      center_(0.0),
      unitOuterNormal_(intersection.centerUnitOuterNormal()),
      area_(geometryHelper.scvfArea(corners_)),
      scvfIndex_(scvfIndex),
      scvIndices_(std:: move(scvIndices)),
      boundary_(boundary)
    , boundaryFlag_{intersection}
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    {
        return center_;
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return area_;
    }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume for spatial param evaluation
    LocalIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    LocalIndexType outsideScvIdx() const
    {
        assert(!boundary());
        return scvIndices_[1];
    }

    //! The local index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    {
        return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_);
    }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    {
        return boundaryFlag_.get();
    }

private:
    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux

#endif
