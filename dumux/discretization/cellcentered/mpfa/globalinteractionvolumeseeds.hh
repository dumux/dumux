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
 * \brief Base class for the global interaction volume seeds of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_GLOBALINTERACTIONVOLUMESEEDS_HH
#define DUMUX_DISCRETIZATION_MPFA_GLOBALINTERACTIONVOLUMESEEDS_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/globalinteractionvolumeseedsbase.hh>
#include "methods.hh"

namespace Dumux
{
//! forward declaration of the actual method-specific implementation
//! By default we simply inherit from the base class
//! Actual implementations for other methods have to be provided below
template<class TypeTag, MpfaMethods method>
class CCMpfaGlobalInteractionVolumeSeedsImplementation : public CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>
{
    using ParentType = CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

public:
    CCMpfaGlobalInteractionVolumeSeedsImplementation(const GridView gridView) : ParentType(gridView)  {}
};

/*!
 * \ingroup Mpfa
 * \brief Base class for the creation and storage of the interaction volume seeds for mpfa methods.
 */
template<class TypeTag>
using CCMpfaGlobalInteractionVolumeSeeds = CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, GET_PROP_VALUE(TypeTag, MpfaMethod)>;

} // end namespace

// the specializations of this class differing from the default have to be included here
#include <dumux/discretization/cellcentered/mpfa/lmethod/globalinteractionvolumeseeds.hh>

#endif
