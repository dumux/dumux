// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMomentumPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH

#include <dumux/discretization/method.hh>

#include "couplingmanager_staggered_cctpfa.hh"

namespace Dumux {

namespace FreeFlowMomentumPorousMediumDetail {

// declaration (specialize for different discretization types)
template<class MDTraits,
         class DiscFFMomentum = typename MDTraits::template SubDomain<0>::GridGeometry::DiscretizationMethod,
         class DiscPM = typename MDTraits::template SubDomain<1>::GridGeometry::DiscretizationMethod
         >
struct FreeFlowMomentumPorousMediumCouplingManagerSelector;

template<class MDTraits>
struct FreeFlowMomentumPorousMediumCouplingManagerSelector<MDTraits, DiscretizationMethods::FCStaggered, DiscretizationMethods::CCTpfa>
{ using type = FFMomentumPMCouplingManagerStaggeredCCTpfa<MDTraits>; };

} // end namespace Detail

template<class MDTraits>
using FreeFlowMomentumPorousMediumCouplingManager = typename FreeFlowMomentumPorousMediumDetail::FreeFlowMomentumPorousMediumCouplingManagerSelector<MDTraits>::type;

} // end namespace Dumux

#endif
