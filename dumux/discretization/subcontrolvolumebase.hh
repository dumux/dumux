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
 * \brief Base class for a sub control volume
 */
#ifndef DUMUX_SUBCONTROLVOLUME_HH
#define DUMUX_SUBCONTROLVOLUME_HH

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Base class for a sub control volume, i.e a part of the control
 *        volume we are making the balance for. Defines the general interface.
 * \tparam Imp the implementation
 * \tparam ScvGeometryTraits traits of this class
 */
template<class Imp, class ScvGeometryTraits>
class SubControlVolumeBase
{
    using Implementation = Imp;
    using GridIndexType = typename ScvGeometryTraits::GridIndexType;
    using LocalIndexType = typename ScvGeometryTraits::LocalIndexType;
    using Scalar = typename ScvGeometryTraits::Scalar;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename ScvGeometryTraits::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = ScvGeometryTraits;

    //! The center of the sub control volume
    GlobalPosition center() const
    {
        return asImp_().center();
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return asImp_().volume();
    }

    //! The index of the dof this scv is embedded in (ccfv)
    GridIndexType dofIndex() const
    {
        return asImp_().dofIndex();
    }

    //! The index of the dof this scv is embedded in (box)
    LocalIndexType localDofIndex() const
    {
        return asImp_().localDofIndex();
    }

    // The position of the dof this scv is embedded in
    GlobalPosition dofPosition() const
    {
        return asImp_().dofPosition();
    }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    {
        return asImp_().elementIndex();
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this);}

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this);}
};

} // end namespace Dumux

#endif
