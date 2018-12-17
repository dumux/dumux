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
 * \ingroup Discretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_SUBCONTROLVOLUMEFACEBASE_HH
#define DUMUX_DISCRETIZATION_SUBCONTROLVOLUMEFACEBASE_HH

#include <utility>
#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Base class for a sub control volume face, i.e a part of the boundary
 *        of a sub control volume we computing a flux on.
 * \tparam Imp the implementation
 * \tparam ScvGeometryTraits traits of this class
 */
template<class Imp, class ScvfGeometryTraits>
class SubControlVolumeFaceBase
{
    using Implementation = Imp;
    using GridIndexType = typename ScvfGeometryTraits::GridIndexType;
    using Scalar = typename ScvfGeometryTraits::Scalar;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename ScvfGeometryTraits::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = ScvfGeometryTraits;

    //! The center of the sub control volume face
    GlobalPosition center() const
    {
        return asImp_().center();
    }

    //! The integration point for flux evaluations in global coordinates
    GlobalPosition ipGlobal() const
    {
        // Return center for now
        return asImp_().ipGlobal();
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return asImp_().area();
    }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return asImp_().boundary();
    }

    //! the unit outward pointing normal on the scv face
    GlobalPosition unitOuterNormal() const
    {
        return asImp_().unitOuterNormal();
    }

    //! index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    {
        return asImp_().insideScvIdx();
    }

    //! index of the outside sub control volume for spatial param evaluation
    //! This results in undefined behaviour if boundary is true
    //! In case of multiple outside scv indices (network grids) an index can be supplied
    GridIndexType outsideScvIdx(int i = 0) const
    {
        return asImp_().outsideScvIdx(i);
    }

    //! The global index of this sub control volume face
    GridIndexType index() const
    {
        return asImp_().index();
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
};

} // end namespace Dumux

#endif
