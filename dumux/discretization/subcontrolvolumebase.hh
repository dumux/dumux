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
#ifndef DUMUX_SUBCONTROLVOLUME_HH
#define DUMUX_SUBCONTROLVOLUME_HH

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief Base class for a sub control volume, i.e a part of the control
 *        volume we are making the balance for. Defines the general interface.
 */
template<class Imp, class G, typename I>
class SubControlVolumeBase
{
    using Implementation = Imp;
    using IndexType = I;
    using Geometry = typename std::decay<G>::type;

    using Scalar = typename Geometry::ctype;
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:

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

    //! The global index of this scv
    IndexType index() const
    {
        return asImp_().index();
    }

    //! The index of the dof this scv is embedded in
    IndexType dofIndex() const
    {
        return asImp_().dofIndex();
    }

    // The position of the dof this scv is embedded in
    GlobalPosition dofPosition() const
    {
        return asImp_().dofPosition();
    }

    //! The global index of the element this scv is embedded in
    IndexType elementIndex() const
    {
        return asImp_().elementIndex();
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this);}

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this);}
};

} // end namespace

#endif
