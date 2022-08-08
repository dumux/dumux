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
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_HH

#include <dumux/discretization/method.hh>

#include "couplingmanager_staggered.hh"
#include "couplingmanager_diamond.hh"

#ifndef DOXYGEN
namespace Dumux::Detail {

template<class Traits>
struct MomentumDiscMethod
{ using type = typename Traits::template SubDomain<0>::GridGeometry::DiscretizationMethod; };

// declaration (specialize for different discretization types)
template<class Traits, class DiscretizationMethod = typename MomentumDiscMethod<Traits>::type>
struct FreeFlowCouplingManagerSelector;

template<class Traits>
struct FreeFlowCouplingManagerSelector<Traits, DiscretizationMethods::FCStaggered>
{ using type = FCStaggeredFreeFlowCouplingManager<Traits>; };

template<class Traits>
struct FreeFlowCouplingManagerSelector<Traits, DiscretizationMethods::FCDiamond>
{ using type = FCDiamondFreeFlowCouplingManager<Traits>; };

} // end namespace Dumux::Detail
#endif


namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for free flow systems
 */
template<class Traits>
using FreeFlowCouplingManager = typename Detail::FreeFlowCouplingManagerSelector<Traits>::type;

} // end namespace Dumux

#endif
