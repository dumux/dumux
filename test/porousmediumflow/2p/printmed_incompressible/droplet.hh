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

    // the contructor
    Droplet(const Scalar initialVolume,
            const Scalar initialRadius,
            const Scalar initialHeight,
            const Scalar initialContactRadius,
            const Scalar initialContactAngle,
            const GlobalPosition initialCenter,
            const std::vector<GridIndexType> dropletElems,
            const std::vector<GridIndexType> dropletDoFs,
            const std::vector<GlobalPosition> dropletDoFPositions)
    : volume_(initialVolume)
    , initialVolume_(initialVolume)
    , radius_(initialRadius)
    , initialRadius_(initialRadius)
    , height_(initialHeight)
    , initialHeight_(initialHeight)
    , initialContactRadius_(initialContactRadius)
    , contactAngle_(initialContactAngle)
    , initialContactAngle_(initialContactAngle)
    , initialCenter_(initialCenter)
    , dropletElems_(dropletElems)
    , dropletDoFs_(dropletDoFs)
    , dropletDoFPositions_(dropletDoFPositions)
    {}


    void update(const Scalar volume,
                const Scalar radius,
                const Scalar height,
                const Scalar contactAngle,)
    {
        volume_ = volume;
        radius_ = radius;
        height_ = height;
        contactAngle_ = contactAngle;
    }

    Scalar radius() const
    {   return radius_; }

    Scalar initialRadius() const
    {   return initialRadius_; }

    Scalar contactAngle() const
    {   return contactAngle_; }

    Scalar initialContactAngle() const
    {   return initialContactAngle_; }

    Scalar volume() const
    {   return volume_; }

    Scalar initialVolume() const
    {   return initialVolume_; }

    Scalar contactRadius() const
    {   return initialContactRadius_; }

    Scalar height() const
    {   return height_; }

    Scalar initialHeight() const
    {   return initialHeight_; }

    GlobalPosition center() const
    {   return initialCenter_; }


    std::vector<GridIndexType> elementIndices() const
    {
        return dropletElems_;
    }

    std::vector<GridIndexType> dofIndices() const
    {
        return dropletDoFs_;
    }

    std::vector<GlobalPosition> dofpositions() const
    {
        return dropletDoFPositions_;
    }

private:
    GlobalPosition initialCenter_ ;
    Scalar radius_ = 0.0;
    Scalar initialRadius_ = 0.0;
    Scalar volume_ = 0.0;
    Scalar initialVolume_ = 0.0;
    Scalar contactAngle_ = 0.0;
    Scalar initialContactAngle_ = 0.0;
    Scalar initialContactRadius_ = 0.0;
    Scalar height_ = 0.0;
    Scalar initialHeight_ = 0.0;
    std::vector<GridIndexType> dropletElems_;
    std::vector<GridIndexType> dropletDoFs_;
    std::vector<GlobalPosition> dropletDoFPositions_;

};

} // end namespace Dumux::PoreNetwork

#endif
