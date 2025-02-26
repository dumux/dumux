// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits to detect periodicity support
 */
#ifndef DUMUX_TYPETRAITS_PERIODIC_HH
#define DUMUX_TYPETRAITS_PERIODIC_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux::Detail {

//! helper struct detecting if a gridGeometry object has a periodicDofMap() function
template<class GG>
using GGPeriodicMapDetector = decltype(std::declval<GG>().periodicDofMap());

template<class GG>
constexpr inline bool hasPeriodicDofMap()
{ return Dune::Std::is_detected<GGPeriodicMapDetector, GG>::value; }

} // end namespace Dumux::Detail

#endif
