// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits to be used with matrix types
 */
#ifndef DUMUX_TYPETRAITS_STATE_HH
#define DUMUX_TYPETRAITS_STATE_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux::Detail {

//! helper struct detecting if a PrimaryVariables object has a state() function
struct hasState
{
    template<class PrimaryVariables>
    auto operator()(PrimaryVariables&& priVars)
    -> decltype(priVars.state())
    {}
};

template<class P>
using DetectPriVarsHaveState = decltype(std::declval<P>().state());

template<class P>
constexpr inline bool priVarsHaveState()
{ return Dune::Std::is_detected<DetectPriVarsHaveState, P>::value; }

} // end namespace Dumux::Detail

#endif
