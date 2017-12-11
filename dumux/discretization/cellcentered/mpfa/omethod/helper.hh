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
 * \brief Helper class for interaction volumes using the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_HELPER_HH

#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag, MpfaMethods m, int dim, int dimWorld>
class MpfaMethodHelper;

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions. This scheme does
 *        not require any additional functionality to be implemented.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/2, /*dimWorld*/2> {};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions embedded in a 3d world.
 *        This scheme does not require any additional functionality to be implemented.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/2, /*dimWorld*/3> {};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in three dimensions.
 *        This scheme does not require any additional functionality to be implemented.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/3, /*dimWorld*/3> {};

} // end namespace

#endif
