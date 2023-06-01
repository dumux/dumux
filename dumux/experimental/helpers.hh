// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Discretization
 * \brief Helpers w.r.t. experimental code.
 */
#ifndef DUMUX_EXPERIMENTAL_HELPERS_HH
#define DUMUX_EXPERIMENTAL_HELPERS_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux::Experimental {

template<class T>
using NewGridVariablesDetector = decltype(std::declval<T>().gridVolVars());

template<class T>
constexpr inline bool hasNewGridVarInterface()
{ return Dune::Std::is_detected<NewGridVariablesDetector, T>::value; }

} // end namespace Dumux::Experimental

#endif
