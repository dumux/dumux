// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \ingroup TimeStepping
 * \brief Concepts in the context of multi-stage time integration.
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_CONCEPTS_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_CONCEPTS_HH

#include <memory>
#include <type_traits>
#include <concepts>

#include <dumux/experimental/new_assembly/dumux/common/variables.hh>

namespace Dumux::Concepts {

//! Concept for multi-stage variables, usable in multi-stage time integration steps.
template<typename T>
concept MultiStageVariables
    = TimeDependentVariables<T>
    and requires(T& t, const T& tConst) {
        typename T::StageVariables;
        typename T::StageParams;
        std::convertible_to<const T&, const typename T::StageVariables&>;

        { t.clearStages() };
        { t.setStageParams(std::shared_ptr<const typename T::StageParams>{}) };
        { tConst.stageParams() } -> std::same_as<const typename T::StageParams&>;
        { tConst.stageVariables(std::size_t{}) } -> std::same_as<const typename T::StageVariables&>;
        { tConst.numStages() } -> std::integral;
};

} // namespace Dumux::Concepts

#endif
