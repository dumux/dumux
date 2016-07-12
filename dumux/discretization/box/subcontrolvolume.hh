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

namespace Dumux
{
template<class G, typename I>
class BoxSubControlVolume : public SubControlVolumeBase<G, I>
{
public:
    // exported types
    using Geometry = G;
    using IndexType = I;

private:
    using Scalar = typename Geometry::ctype;
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    // the contructor in the box case
    BoxSubControlVolume(Geometry&& geometry,
                        IndexType scvIdx,
                        IndexType elementIndex,
                        IndexType dofIndex)
    : SubControlVolumeBase<G, I>(std::move(geometry), std::move(elementIndex)),
      scvIdx_(std::move(scvIdx)),
      dofIndex_(std::move(dofIndex)) {}

    //! The local index of this scv
    IndexType index() const
    {
        return scvIdx_;
    }

    IndexType dofIndex() const
    {
        return dofIndex_;
    }

    GlobalPosition dofPosition() const
    {
        return this->geometry().corner(scvIdx_);
    }

private:
    IndexType scvIdx_;
    IndexType dofIndex_;
};

} // end namespace

#endif
