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

#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/common/math.hh>

namespace Dumux
{
template<class G, typename I>
class BoxSubControlVolume : public SubControlVolumeBase<BoxSubControlVolume<G, I>, G, I>
{
    using ParentType = SubControlVolumeBase<BoxSubControlVolume<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    enum { dim = Geometry::mydimension };
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    // the contructor in the box case
    BoxSubControlVolume(std::vector<GlobalPosition>&& corners,
                        Scalar volume,
                        IndexType scvIdx,
                        IndexType elementIndex,
                        IndexType dofIndex)
    : ParentType(),
      corners_(std::move(corners)),
      center_(0.0),
      volume_(volume),
      elementIndex_(std::move(elementIndex)),
      scvIdx_(std::move(scvIdx)),
      dofIndex_(std::move(dofIndex))
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
        return Geometry(Dune::GeometryType(Dune::GeometryType::cube, dim), corners_);
    }

    //! The global index of this scv
    IndexType index() const
    {
        return scvIdx_;
    }

    //! The index of the dof this scv is embedded in
    IndexType dofIndex() const
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
    IndexType elementIndex() const
    {
        return elementIndex_;
    }

    //! Return the corner for the given local index
    GlobalPosition corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

private:
    std::vector<GlobalPosition> corners_;
    GlobalPosition center_;
    Scalar volume_;
    IndexType elementIndex_;
    IndexType scvIdx_;
    IndexType dofIndex_;
};

} // end namespace

#endif
