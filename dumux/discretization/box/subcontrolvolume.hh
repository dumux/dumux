// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for a sub control volume
 */
#ifndef DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH

#include <dune/common/version.hh>

#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/common/math.hh>

namespace Dumux
{
template<class ScvGeometryTraits>
class BoxSubControlVolume : public SubControlVolumeBase<BoxSubControlVolume<ScvGeometryTraits>, ScvGeometryTraits>
{
    using ParentType = SubControlVolumeBase<BoxSubControlVolume<ScvGeometryTraits>, ScvGeometryTraits>;
    using Geometry = typename ScvGeometryTraits::Geometry;
    using GridIndexType = typename ScvGeometryTraits::GridIndexType;
    using LocalIndexType = typename ScvGeometryTraits::LocalIndexType;
    using Scalar = typename ScvGeometryTraits::Scalar;
    using GlobalPosition = typename ScvGeometryTraits::GlobalPosition;
    using CornerStorage = typename ScvGeometryTraits::CornerStorage;
    enum { dim = Geometry::mydimension };

public:
    //! state the traits public and thus export all types
    using Traits = ScvGeometryTraits;

    //! The default constructor
    BoxSubControlVolume() = default;

    // the contructor in the box case
    template<class GeometryHelper>
    BoxSubControlVolume(const GeometryHelper& geometryHelper,
                        LocalIndexType scvIdx,
                        GridIndexType elementIndex,
                        GridIndexType dofIndex)
    : corners_(geometryHelper.getScvCorners(scvIdx)),
      center_(0.0),
      volume_(geometryHelper.scvVolume(corners_)),
      elementIndex_(elementIndex),
      scvIdx_(scvIdx),
      dofIndex_(dofIndex)
    {
        // compute center point
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! The center of the sub control volume
    GlobalPosition center() const
    {
        return center_;
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return volume_;
    }

    //! The geometry of the sub control volume
    // e.g. for integration
    Geometry geometry() const
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Geometry(Dune::GeometryTypes::cube(dim), corners_);
#else
        return Geometry(Dune::GeometryType(Dune::GeometryType::cube, dim), corners_);
#endif
    }

    //! The global index of this scv
    LocalIndexType indexInElement() const
    {
        return scvIdx_;
    }

    //! The index of the dof this scv is embedded in
    GridIndexType dofIndex() const
    {
        return dofIndex_;
    }

    // The position of the dof this scv is embedded in
    GlobalPosition dofPosition() const
    {
        // The corner list is defined such that the first entry is the vertex itself
        return corner(0);
    }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    {
        return elementIndex_;
    }

    //! Return the corner for the given local index
    GlobalPosition corner(LocalIndexType localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

private:
    CornerStorage corners_;
    GlobalPosition center_;
    Scalar volume_;
    GridIndexType elementIndex_;
    LocalIndexType scvIdx_;
    GridIndexType dofIndex_;
};

} // end namespace

#endif
