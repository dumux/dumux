// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_HH

#include <dumux/discretization/method.hh>

#include "couplingmanager_base.hh"
#include "couplingmanager_staggered_cctpfa.hh"
namespace Dumux {

namespace FreeFlowPorousMediumDetail {

// declaration (specialize for different discretization types)
template<class MDTraits,
         class DiscFFMomentum = typename MDTraits::template SubDomain<FreeFlowPorousMediumDetail::freeFlowMomentumIndex>::GridGeometry::DiscretizationMethod,
         class DiscFFMass = typename MDTraits::template SubDomain<FreeFlowPorousMediumDetail::freeFlowMassIndex>::GridGeometry::DiscretizationMethod,
         class DiscPM = typename MDTraits::template SubDomain<FreeFlowPorousMediumDetail::porousMediumIndex>::GridGeometry::DiscretizationMethod
         >
struct FreeFlowPorousMediumCouplingManagerSelector;

template<class MDTraits>
struct FreeFlowPorousMediumCouplingManagerSelector<MDTraits, DiscretizationMethods::FCStaggered, DiscretizationMethods::CCTpfa, DiscretizationMethods::CCTpfa>
{ using type = FreeFlowPorousMediumCouplingManagerStaggeredCCTpfa<MDTraits>; };

} // end namespace Detail

template<class MDTraits>
using FreeFlowPorousMediumCouplingManager = typename FreeFlowPorousMediumDetail::FreeFlowPorousMediumCouplingManagerSelector<MDTraits>::type;

} // end namespace Dumux

#endif
