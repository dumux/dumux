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
 * \ingroup CCMpfaDiscretization
 * \brief A helper class to fill the flux variable caches used in the flux constitutive laws
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/porousmediumflow/fluxvariablescachefiller.hh>
#warning "This header is deprecated and will be removed after 3.1, use dumux/porousmediumflow/fluxvariablescachefiller.hh"

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to fill the flux variables caches within
 *        the interaction volume around a given sub-control volume face.
 */
template<class TypeTag>
using CCMpfaFluxVariablesCacheFiller [[deprecated("This class has been renamed to PorousMediumFluxVariablesCacheFiller and will be removed after 3.1")]]
= PorousMediumFluxVariablesCacheFiller<TypeTag>;

} // end namespace Dumux

#endif
