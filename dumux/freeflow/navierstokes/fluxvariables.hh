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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/freeflow/navierstokes/staggered/fluxvariables.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableInertiaTerms);
}

// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class NavierStokesFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using NavierStokesFluxVariables = NavierStokesFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;


} // end namespace

#endif
