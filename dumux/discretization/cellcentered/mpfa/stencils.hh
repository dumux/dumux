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
 * \brief Implements the notion of stencils for cell-centered models
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_STENCILS_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_STENCILS_HH

#include <dumux/discretization/cellcentered/stencils.hh>

namespace Dumux
{
//! Method-specific implementation of the stencils vector
//! By default, we inherit from the symmetric stencils vector class
//! For non-symmetric mpfa methods a corresponding specialization has to be provided below
template<class TypeTag, MpfaMethods method>
class CCMpfaStencilsVectorImplementation : public CCSymmetricStencilsVector<TypeTag> {};

template<class TypeTag>
using CCMpfaStencilsVector = CCMpfaStencilsVectorImplementation<TypeTag, GET_PROP_VALUE(TypeTag, MpfaMethod)>;

} // end namespace

#endif
