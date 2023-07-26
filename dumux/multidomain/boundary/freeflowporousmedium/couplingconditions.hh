// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling conditions specialized for different discretization schemes
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_HH

#include <dumux/discretization/method.hh>

#include "couplingconditions_staggered_cctpfa.hh"

namespace Dumux {

#ifndef DOXYGEN
namespace FreeFlowPorousMediumDetail {

// declaration (specialize for different discretization types)
template<class MDTraits, class CouplingManager,
         class DiscFFMomentum = typename MDTraits::template SubDomain<CouplingManager::freeFlowMomentumIndex>::GridGeometry::DiscretizationMethod,
         class DiscFFMass = typename MDTraits::template SubDomain<CouplingManager::freeFlowMassIndex>::GridGeometry::DiscretizationMethod,
         class DiscPM = typename MDTraits::template SubDomain<CouplingManager::porousMediumIndex>::GridGeometry::DiscretizationMethod
         >
struct FreeFlowPorousMediumCouplingConditionsSelector;

template<class MDTraits, class CouplingManager>
struct FreeFlowPorousMediumCouplingConditionsSelector<MDTraits, CouplingManager, DiscretizationMethods::FCStaggered, DiscretizationMethods::CCTpfa, DiscretizationMethods::CCTpfa>
{ using type = FFPMCouplingConditionsStaggeredCCTpfa<MDTraits, CouplingManager>; };

} // end namespace FreeFlowPorousMediumDetail
#endif // DOXYGEN

template<class MDTraits, class CouplingManager>
using FreeFlowPorousMediumCouplingConditions = typename FreeFlowPorousMediumDetail::FreeFlowPorousMediumCouplingConditionsSelector<MDTraits,CouplingManager>::type;

} // end namespace Dumux

#endif
