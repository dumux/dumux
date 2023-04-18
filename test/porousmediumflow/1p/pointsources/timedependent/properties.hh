// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
