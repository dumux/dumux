// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits for problem classes
 */
#ifndef DUMUX_TYPETRAITS_PROBLEM_HH
#define DUMUX_TYPETRAITS_PROBLEM_HH

#include <type_traits>
#include <dumux/discretization/method.hh>

namespace Dumux {

// forward declare
namespace Detail {
template<class Problem, class DiscretizationMethod>
struct ProblemTraits;

template<class Problem>
struct ProblemGridGeometryHelper
{
private:
    template<class P>
    static auto deduce(int) -> std::decay_t<decltype(std::declval<const P&>().gridGeometry())>;

    template<class P>
    static auto deduce(long) -> std::decay_t<decltype(std::declval<const P&>().gridDiscretization())>;

public:
    using type = decltype(deduce<Problem>(0));
};

template<class Problem>
using ProblemGridGeometry = typename ProblemGridGeometryHelper<Problem>::type;
} // end namespace Detail

/*!
 * \ingroup Typetraits
 * \brief Type traits for problem classes.
 */
template<class Problem>
struct ProblemTraits
{
    using GridGeometry = Detail::ProblemGridGeometry<Problem>;
    using BoundaryTypes = typename Detail::template ProblemTraits<Problem, typename GridGeometry::DiscretizationMethod>::BoundaryTypes;
};

} // end namespace Dumux

#endif
