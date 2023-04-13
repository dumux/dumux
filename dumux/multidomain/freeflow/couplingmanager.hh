// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Freeflow coupling managers (Navier-Stokes mass-momentum coupling)
 */
#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_COUPLING_MANAGER_HH

#include <dumux/discretization/method.hh>

#include "typetraits.hh"
#include "couplingmanager_staggered.hh"
#include "couplingmanager_cvfe.hh"

#ifndef DOXYGEN
namespace Dumux::Detail {

// declaration (specialize for different discretization types)
template<class Traits, class DiscretizationMethod = typename MomentumDiscretizationMethod<Traits>::type>
struct FreeFlowCouplingManagerSelector;

template<class Traits>
struct FreeFlowCouplingManagerSelector<Traits, DiscretizationMethods::FCStaggered>
{ using type = FCStaggeredFreeFlowCouplingManager<Traits>; };

template<class Traits, class D>
struct FreeFlowCouplingManagerSelector<Traits, DiscretizationMethods::CVFE<D>>
{ using type = CVFEFreeFlowCouplingManager<Traits>; };

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
