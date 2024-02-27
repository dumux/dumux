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
#ifndef DUMUX_DROPLET_HH
#define DUMUX_DROPLET_HH

#include <dune/geometry/affinegeometry.hh>
#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>

namespace Dumux{

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Default traits class
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct DropDefaultScvGeometryTraits
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
         class T = DropDefaultScvGeometryTraits<GV> >
class Droplet
{
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    Droplet() = default;

    // the contructor in the box case
    void Droplet(const Scalar volume,
                const Scalar radius,
                const Scalar height,
                const Scalar contactAngle,)
    {
        intrinsicVolume_ = volume;
        intrinsicRadius_ = radius;
        height_ = height;
        contactAngle_ = contactAngle;
    }


    update(const GridIndexType dofIndex,
            const LocalIndexType scvLocalIdx,
            const GridIndexType elementIndex,
            const GlobalPosition& center,
            Scalar contactRadius)
    : dropDofIndex_(dofIndex)
    , dropScvLocalIndex_(scvLocalIdx)
    , dropElementIndex_(elementIndex)
    , center_(center)
    , initialCenter_(center)
    , stickRadius_(contactRadius)
    , volume_(0.0)
    , intrinsicVolume_(0.0)
    , radius_(0.0)
    , intrinsicRadius_(0.0)
    , height_(0.0)
    , contactAngle_(0.0)
    {}


GridIndexType dofIndex() const
    {
        return dropDofIndex_;
    }

    Scalar radius() const
    {
        if (!isMerged_)
            return intrinsicRadius_;

        return radius_;
    }

    Scalar intrinsicRadius() const
    {
         return intrinsicRadius_;
    }

    Scalar contactAngle() const
    {
        return contactAngle_;
    }

    Scalar volume() const
    {
        if (!isMerged_)
            return intrinsicVolume_;

         return volume_;
    }

    Scalar intrinsicVolume() const
    {
         return intrinsicVolume_;
    }

    Scalar stickRadius() const
    {
        return stickRadius_;
    }

    Scalar initialStickRadius() const
    {
        return poreRadius_;
    }

    Scalar height() const
    {
        return height_;
    }

    GlobalPosition center() const
    {
        return center_;
    }
    GlobalPosition initialCenter() const
    {
        return initialCenter_;
    }

    bool isOnBoundary() const
    {
        return isOnOutletBoundary_;
    }

    GridIndexType elementIdx() const
    {
        return dropElementIndex_;
    }

    GridIndexType LocalDofIdx() const
    {
        return dropScvLocalIndex_;
    }

private:
    CornerStorage corners_;
    GridIndexType dropDofIndex_;
    GridIndexType dropElementIndex_;
    GridIndexType dropScvLocalIndex_;
    GlobalPosition initialCenter_;
    GlobalPosition center_ ;
    Scalar radius_ = 0.0;
    Scalar intrinsicRadius_ = 0.0;
    Scalar volume_ = 0.0;
    Scalar intrinsicVolume_ = 0.0;
    Scalar contactAngle_ = 0.0;
    Scalar stickRadius_ = 0.0;
    Scalar height_ = 0.0;
    bool isOnOutletBoundary_;
};

} // end namespace Dumux::PoreNetwork

#endif
