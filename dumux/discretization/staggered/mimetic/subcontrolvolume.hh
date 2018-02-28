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
#ifndef DUMUX_DISCRETIZATION_MIMETIC_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_MIMETIC_SUBCONTROLVOLUME_HH

#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/common/optional.hh>

namespace Dumux
{
template<class G, typename I>
class MimeticSubControlVolume : public SubControlVolumeBase<MimeticSubControlVolume<G, I>, G, I>
{
    using ParentType = SubControlVolumeBase<MimeticSubControlVolume<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    // the default constructor
    MimeticSubControlVolume() = default;

    // the contructor in the cc case
    MimeticSubControlVolume(Geometry&& geometry,
                       IndexType elementIndex)
    : ParentType(), geometry_(std::forward<Geometry>(geometry)), elementIndex_(elementIndex),
                    center_(geometry.center()), volume_(geometry.volume())
    {
    }

    //! The copy constrcutor
    MimeticSubControlVolume(const MimeticSubControlVolume& other) = default;

    //! The move constrcutor
    MimeticSubControlVolume(MimeticSubControlVolume&& other) = default;

    //! The copy assignment operator
    MimeticSubControlVolume& operator=(const MimeticSubControlVolume& other)
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(other.geometry_.value());
        elementIndex_ = other.elementIndex_;
        center_ = geometry().center();
        volume_ = geometry().volume();
        return *this;
    }

    //! The move assignment operator
    MimeticSubControlVolume& operator=(MimeticSubControlVolume&& other)
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(other.geometry_.value());
        elementIndex_ = std::move(other.elementIndex_);
        center_ = geometry().center();
        volume_ = geometry().volume();
        return *this;
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
    const Geometry& geometry() const
    {
        assert((geometry_));
        return geometry_.value();
    }

    //! The index of the dof this scv is embedded in (the global index of this scv)
    IndexType dofIndex() const
    {
        return elementIndex();
    }

    //! The global index of this scv
    IndexType indexInElement() const
    {
        return 0;
    }

    // The position of the dof this scv is embedded in
    GlobalPosition dofPosition() const
    {
        return center_;
    }

    //! The global index of the element this scv is embedded in
    IndexType elementIndex() const
    {
        return elementIndex_;
    }

    //! Return the corner for the given local index
    GlobalPosition corner(unsigned int localIdx) const
    {
        assert(localIdx < geometry().corners() && "provided index exceeds the number of corners");
        return geometry().corner(localIdx);
    }

    void setCellCenter(GlobalPosition center)
    {
        center_ = center;
    }

    void setCellVolume(Scalar volume)
    {
        volume_ = volume;
    }

private:
    // Work around the fact that geometry is not default constructible
    Optional<Geometry> geometry_;
    IndexType elementIndex_;
    GlobalPosition center_;
    Scalar volume_;
};
} // end namespace

#endif
