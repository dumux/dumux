// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
