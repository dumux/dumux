// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Traits shared by the multi-domain assemblers
 */
#ifndef DUMUX_MULTIDOMAIN_ASSEMBLY_TRAITS_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLY_TRAITS_HH

#include <type_traits>
#include <utility>

#include <dune/common/std/type_traits.hh>

namespace Dumux {

namespace Detail {

//! helper struct detecting if sub-problem has a constraints() function
template<class P>
using SubProblemConstraintsDetector = decltype(std::declval<P>().constraints());

template<class P>
constexpr inline bool hasSubProblemGlobalConstraints()
{ return Dune::Std::is_detected<SubProblemConstraintsDetector, P>::value; }

} // end namespace Detail

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Type trait that is specialized for coupling manager supporting multithreaded assembly
 * \note A coupling manager implementation that wants to enable multithreaded assembly has to specialize this trait
 */
template<class CM>
struct CouplingManagerSupportsMultithreadedAssembly : public std::false_type
{};

} // end namespace Dumux

#endif
