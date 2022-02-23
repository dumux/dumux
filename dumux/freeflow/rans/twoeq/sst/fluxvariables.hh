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
 * \ingroup SSTModel
 * \copydoc Dumux::SSTFluxVariables
 */
#ifndef DUMUX_SST_FLUXVARIABLES_HH
#define DUMUX_SST_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/rans/twoeq/sst/staggered/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseFluxVariables, class DiscretizationMethod>
class SSTFluxVariablesImpl;

/*!
 * \ingroup SSTModel
 * \brief The flux variables class for the SST model.
          This is a convenience alias for that actual,
          discretization-specific flux variables.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag, class BaseFluxVariables>
using SSTFluxVariables = SSTFluxVariablesImpl<TypeTag, BaseFluxVariables, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace

#endif
