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
 * \ingroup PoreNetworkDiscretization
 * \brief the sub control volume for pore networks
 */
#ifndef DUMUX_DISCRETIZATION_PNM_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_PNM_SUBCONTROLVOLUME_HH

#include <dune/geometry/affinegeometry.hh>
#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Default traits class
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PNMDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CornerStorage = std::array<GlobalPosition, 2>;
    using Geometry = Dune::AffineGeometry<Scalar, 1, dimWorld>;
};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief the sub control volume for porenetworks
 * \tparam GV the type of the grid view
 * \tparam T the scv geometry traits
 */
template<class GV,
         class T = PNMDefaultScvGeometryTraits<GV> >
class PNMSubControlVolume
: public Dumux::SubControlVolumeBase<PNMSubControlVolume<GV, T>, T>
{
    using ThisType = PNMSubControlVolume<GV, T>;
    using ParentType = Dumux::SubControlVolumeBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    PNMSubControlVolume() = default;

    // the contructor in the box case
    template<class Corners>
    PNMSubControlVolume(GridIndexType dofIndex,
                        LocalIndexType scvIdx,
                        GridIndexType elementIndex,
                        Corners&& corners,
                        const Scalar volume)
    : center_((corners[0]+corners[1])/2.0),
      corners_(std::forward<CornerStorage>(corners)),
      volume_(volume),
      elementIndex_(elementIndex),
      localDofIdx_(scvIdx),
      dofIndex_(dofIndex)
    {}

    //! The center of the sub control volume (return pore center).
    //! Be aware that this is not the pore-body center! Use dofPosition() for the latter!
    const GlobalPosition& center() const
    { return center_; }

    //! The volume of the sub control volume (part of a pore)
    Scalar volume() const
    { return volume_; }

    //! The element-local index of the dof this scv is embedded in
    LocalIndexType localDofIndex() const
    { return localDofIdx_; }

    //! The element-local index of this scv.
    LocalIndexType indexInElement() const
    { return localDofIdx_; }

    //! The index of the dof this scv is embedded in
    GridIndexType dofIndex() const
    { return dofIndex_; }

    // The position of the dof this scv is embedded in
    const GlobalPosition& dofPosition() const
    { return corners_[0]; }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    { return elementIndex_; }

    //! Return the corner for the given local index
    const GlobalPosition& corner(LocalIndexType localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume e.g. for integration
    Geometry geometry() const
    {
        return Geometry(Dune::GeometryTypes::simplex(1), corners_);
    }

private:
    GlobalPosition center_;
    CornerStorage corners_;
    Scalar volume_;
    GridIndexType elementIndex_;
    LocalIndexType localDofIdx_;
    GridIndexType dofIndex_;
};

} // end namespace Dumux::PoreNetwork

#endif
