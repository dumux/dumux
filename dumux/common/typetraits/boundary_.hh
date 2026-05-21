// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits for boundary handling
 */
#ifndef DUMUX_TYPETRAITS_BOUNDARY_HH
#define DUMUX_TYPETRAITS_BOUNDARY_HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>

namespace Dumux::Detail {

//! helper struct detecting if a problem has boundaryFlux function
template<class P, class FVG, class EV, class EFVC, class IP>
using BoundaryFluxFunctionDetector = decltype(
    std::declval<P>().boundaryFlux(std::declval<FVG>(), std::declval<EV>(), std::declval<EFVC>(), std::declval<IP>())
);

template<class P, class FVG, class EV, class EFVC, class IP>
constexpr inline bool hasProblemBoundaryFluxFunction()
{ return Dune::Std::is_detected<BoundaryFluxFunctionDetector, P, FVG, EV, EFVC, IP>::value; }

//! helper struct detecting if a problem has boundaryTypes function
template<class P, class FVG, class FACE>
using BoundaryTypesForFaceFunctionDetector = decltype(
    std::declval<P>().boundaryTypes(std::declval<FVG>(), std::declval<FACE>())
);

template<class P, class FVG>
constexpr inline bool hasProblemBoundaryTypesForFaceFunction()
{ return Dune::Std::is_detected<BoundaryTypesForFaceFunctionDetector, P, FVG, typename FVG::BoundaryFace>::value; }

//! helper struct to determine BoundaryTypes related to provided boundaryTypes(...) function of the problem
template<class P, class FVG, bool newInterface = hasProblemBoundaryTypesForFaceFunction<P, FVG>()>
struct BoundaryTypes;

template<class P, class FVG>
struct BoundaryTypes<P, FVG, false>
{ using type = std::decay_t< decltype( std::declval<P>().boundaryTypes(std::declval<typename FVG::Element>(), std::declval<typename FVG::SubControlVolume>()) ) >; };

//! helper struct to determine boundary type
template<class P, class FVG>
struct BoundaryTypes<P, FVG, true>
{ using type = std::decay_t< decltype( std::declval<P>().boundaryTypes(std::declval<FVG>(), std::declval<typename FVG::BoundaryFace>()) ) >; };

} // end namespace Dumux::Detail

#endif
