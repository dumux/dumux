// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMassPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMASSPM_COUPLINGMANAGER_HH

#include <dumux/discretization/method.hh>

#include "couplingmanager_staggered_cctpfa.hh"

namespace Dumux {

namespace FreeFlowMassPorousMediumDetail {

// declaration (specialize for different discretization types)
template<class MDTraits,
         class DiscFFMass = typename MDTraits::template SubDomain<0>::GridGeometry::DiscretizationMethod,
         class DiscPM = typename MDTraits::template SubDomain<1>::GridGeometry::DiscretizationMethod
         >
struct FreeFlowMassPorousMediumCouplingManagerSelector;

template<class MDTraits>
struct FreeFlowMassPorousMediumCouplingManagerSelector<MDTraits, DiscretizationMethods::CCTpfa, DiscretizationMethods::CCTpfa>
{ using type = FFMassPMCouplingManagerStaggeredCCTpfa<MDTraits>; };

} // end namespace Detail

template<class MDTraits>
using FreeFlowMassPorousMediumCouplingManager = typename FreeFlowMassPorousMediumDetail::FreeFlowMassPorousMediumCouplingManagerSelector<MDTraits>::type;

} // end namespace Dumux

#endif
