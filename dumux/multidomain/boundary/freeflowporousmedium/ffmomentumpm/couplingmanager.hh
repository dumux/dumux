// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling managers specialized for different discretization schemes for momentum coupling
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FFPM_FFMMOMENTUMPM_COUPLINGMANAGER_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

#include "couplingmanager_staggered_cctpfa.hh"

namespace Dumux {

#ifndef DOXYGEN
namespace FreeFlowMomentumPorousMediumDetail {

// declaration (specialize for different discretization types)
template<class MDTraits,
         class DiscFFMomentum = typename GetPropType<typename MDTraits::template SubDomainTypeTag<0>, Properties::GridGeometry>::DiscretizationMethod,
         class DiscPM = typename GetPropType<typename MDTraits::template SubDomainTypeTag<1>, Properties::GridGeometry>::DiscretizationMethod
         >
struct FreeFlowMomentumPorousMediumCouplingManagerSelector;

template<class MDTraits>
struct FreeFlowMomentumPorousMediumCouplingManagerSelector<MDTraits, DiscretizationMethods::FCStaggered, DiscretizationMethods::CCTpfa>
{ using type = FFMomentumPMCouplingManagerStaggeredCCTpfa<MDTraits>; };

} // end namespace FreeFlowMomentumPorousMediumDetail
#endif // DOXYGEN

template<class MDTraits>
using FreeFlowMomentumPorousMediumCouplingManager = typename FreeFlowMomentumPorousMediumDetail::FreeFlowMomentumPorousMediumCouplingManagerSelector<MDTraits>::type;

} // end namespace Dumux

#endif
