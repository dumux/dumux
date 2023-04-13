// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test problem for the one-phase model:
 * Water is injected in one single point in the middle of the domain.
 */

#ifndef DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_PROPERTIES_HH
#define DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_PROPERTIES_HH

#include "../timeindependent/properties.hh"
#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePSingularityTimeDependentCCTpfa { using InheritsFrom = std::tuple<OnePSingularityCCTpfa>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSingularityTimeDependentCCTpfa> { using type = OnePSingularityProblemTimeDependent<TypeTag>; };

// point source
template<class TypeTag>
struct PointSource<TypeTag, TTag::OnePSingularityTimeDependentCCTpfa> { using type = SolDependentPointSource<TypeTag>; };
} // end namespace Dumux::Properties

#endif
