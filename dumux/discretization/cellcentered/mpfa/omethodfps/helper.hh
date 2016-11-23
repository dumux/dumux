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
 * \brief Helper class to get the required information on an interaction volume.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFAO_FPS_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFAO_FPS_HELPER_HH

#include <dumux/discretization/cellcentered/mpfa/omethod/helper.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions.
 */
template<class TypeTag, int dim, int dimWorld>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethodFps, dim, dimWorld> : public MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, dim, dimWorld>
{};

} // end namespace

#endif
