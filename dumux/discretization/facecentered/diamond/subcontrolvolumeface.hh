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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FreeFlowDiamondSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUMEFACE_HH

#include <array>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>


#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct FaceCenteredDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::SimplexTopology< mydim >::type::id;
        };
    };

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar>>;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;

    static Dune::GeometryType geometryType()
    {
        return Dune::GeometryTypes::simplex(dim-1);
    }
};

template<class GridView, class T = FaceCenteredDefaultScvfGeometryTraits<GridView>>
class FaceCenteredDiamondSubControlVolumeFace
{
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    using GlobalPosition = typename T::GlobalPosition;

    FaceCenteredDiamondSubControlVolumeFace() = default;

    FaceCenteredDiamondSubControlVolumeFace(const CornerStorage& corners,
                                            const GlobalPosition& normal,
                                            const std::array<LocalIndexType, 2> localScvIndices,
                                            const std::array<GridIndexType, 2> globalScvIndices,
                                            const LocalIndexType localScvfIdx,
                                            const GridIndexType globalScvfIdx,
                                            const bool boundary = false)
    : corners_(corners)
    , geometry_(Geometry(T::geometryType(), corners_))
    , center_(geometry_.value().center())
    , ipGlobal_(geometry_.value().center())
    , unitOuterNormal_(normal)
    , localScvIndices_(localScvIndices)
    , globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , area_(geometry_.value().volume())
    , globalScvfIdx_(globalScvfIdx)
    , boundary_(boundary) {}

    template<class Intersection>
    FaceCenteredDiamondSubControlVolumeFace(const CornerStorage& corners,
                                            const Intersection& intersection,
                                            const std::array<LocalIndexType, 2> localScvIndices,
                                            const std::array<GridIndexType, 2> globalScvIndices,
                                            const LocalIndexType localScvfIdx,
                                            const GridIndexType globalScvfIdx,
                                            const bool boundary = true)
    : corners_(corners)
    , geometry_(Geometry(intersection.geometry().type(), corners_))
    , center_(geometry_.value().center())
    , ipGlobal_(geometry_.value().center())
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , localScvIndices_(localScvIndices)
    , globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , area_(geometry_.value().volume())
    , globalScvfIdx_(globalScvfIdx)
    , boundary_(boundary) {}

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point of the sub control volume face
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! Index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    { return globalScvIndices_[0]; }

    //! index of the outside sub control volume for spatial param evaluation
    GridIndexType outsideScvIdx() const
    { return globalScvIndices_[1]; }

    GridIndexType index() const
    { return globalScvfIdx_; }

    LocalIndexType localIndex() const
    { return localScvfIdx_; }

    bool boundary() const
    { return boundary_; }

    Scalar area() const
    { return area_; }

    std::size_t scvfIdxWithCommonEntity() const
    { return scvfIdxWithCommonEntity_; }

private:
    CornerStorage corners_;
    std::optional<Geometry> geometry_;
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    std::array<LocalIndexType, 2> localScvIndices_;
    std::array<GridIndexType, 2> globalScvIndices_;
    LocalIndexType localScvfIdx_;
    Scalar area_;
    GridIndexType globalScvfIdx_;
    GridIndexType scvfIdxWithCommonEntity_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH
