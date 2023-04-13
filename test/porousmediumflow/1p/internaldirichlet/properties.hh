// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test for internal Dirichlet constraints
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_PROPERTIES_HH

#include <test/porousmediumflow/1p/incompressible/properties.hh>
#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePInternalDirichlet {};
struct OnePInternalDirichletTpfa { using InheritsFrom = std::tuple<OnePInternalDirichlet, OnePIncompressibleTpfa>; };
struct OnePInternalDirichletBox { using InheritsFrom = std::tuple<OnePInternalDirichlet, OnePIncompressibleBox>; };
} // end namespace TTag

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePInternalDirichlet>
{ using type = OnePTestProblemInternalDirichlet<TypeTag>; };

} // end namespace Dumux::Properties

#endif
