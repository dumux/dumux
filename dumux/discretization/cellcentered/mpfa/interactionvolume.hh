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
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUME_HH

#include "methods.hh"

namespace Dumux
{
// forward declaration of the base class
template<class TypeTag, MpfaMethods MpfaMethod>
class CCMpfaInteractionVolumeImplementation;

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa method
 */
template<class TypeTag>
using CCMpfaInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, GET_PROP_VALUE(TypeTag, MpfaMethod)>;

} // end namespace

// the specializations of this class for the available methods have to be included here
#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>

#endif
