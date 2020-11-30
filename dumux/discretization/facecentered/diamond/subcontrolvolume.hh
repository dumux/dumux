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
 * \copydoc Dumux::FaceCenteredDiamondSubControlVolume
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUME_HH

#include <array>
#include <utility>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/indextraits.hh>


#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct FaceCenteredDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    template <class ct>
    using ScvMLGTraits = Dune::MultiLinearGeometryTraits<ct>;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, ScvMLGTraits<Scalar>>;

    using CornerStorage = typename ScvMLGTraits<Scalar>::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;

    static Dune::GeometryType geometryType()
    {
        if constexpr (dim == 3)
            return Dune::GeometryTypes::pyramid;
        else
            return Dune::GeometryTypes::simplex(dim);
    }
};


template<class GridView, class T = FaceCenteredDefaultScvGeometryTraits<GridView>>
class FaceCenteredDiamondSubControlVolume
{
    using GlobalPosition = typename T::GlobalPosition;
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    FaceCenteredDiamondSubControlVolume() = default;

    FaceCenteredDiamondSubControlVolume(const CornerStorage& corners,
                                        const GlobalPosition& dofPosition,
                                        const GridIndexType globalIndex,
                                        const LocalIndexType indexInElement,
                                        const GridIndexType dofIdx,
                                        const GridIndexType eIdx,
                                        const bool boundary)
    : corners_(corners)
    , geometry_(Geometry(T::geometryType(), corners_))
    , center_(geometry_.value().center())
    , dofPosition_(dofPosition)
    , volume_(geometry_.value().volume())
    , globalIndex_(globalIndex)
    , indexInElement_(indexInElement)
    , dofIdx_(dofIdx)
    , eIdx_(eIdx)
    , boundary_(boundary)
    {}

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The position of the degree of freedom
    const GlobalPosition& dofPosition() const
    { return dofPosition_; }

    Scalar volume() const
    { return volume_; }

    GridIndexType dofIndex() const
    { return dofIdx_; }

    GridIndexType index() const
    { return globalIndex_; }

    GridIndexType elementIndex() const
    { return eIdx_; }

    LocalIndexType indexInElement() const
    { return indexInElement_; }

    LocalIndexType localDofIndex() const
    { return indexInElement_; }

    bool boundary() const
    { return boundary_; }

private:
    CornerStorage corners_;
    std::optional<Geometry> geometry_;
    GlobalPosition center_;
    GlobalPosition dofPosition_;
    Scalar volume_;
    GridIndexType globalIndex_;
    LocalIndexType indexInElement_;
    GridIndexType dofIdx_;
    GridIndexType eIdx_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUME_HH
