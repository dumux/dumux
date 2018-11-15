// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup RichardsModel
 * \brief A Richards model newton solver.
 */
#ifndef DUMUX_RICHARDS_PRIVAR_SWITCH_NEWTON_SOLVER_HH
#define DUMUX_RICHARDS_PRIVAR_SWITCH_NEWTON_SOLVER_HH

#include <dumux/common/properties.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/nonlinear/privarswitchnewtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

template <class TypeTag, class Assembler, class LinearSolver, bool enableWaterDiffusionInAir>
class RichardsPrivarSwitchNewtonSolverImplementation;
/*!
 * \ingroup RichardsModel
 * \brief A base for the richards newton solver which derives from the right base newton solver.
  */
template <class TypeTag, class Assembler, class LinearSolver>
using RichardsPrivarSwitchNewtonSolver = RichardsPrivarSwitchNewtonSolverImplementation <TypeTag, Assembler, LinearSolver, GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir)>;

/*!
 * \ingroup RichardsModel
 * \brief the case without a primary variables switch
  */
template <class TypeTag, class Assembler, class LinearSolver>
class RichardsPrivarSwitchNewtonSolverImplementation<TypeTag, Assembler, LinearSolver, false> : public NewtonSolver<Assembler, LinearSolver>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver>;
public:
    using ParentType::ParentType;
};
/*!
 * \ingroup RichardsModel
 * \brief the case with switchable primary variables
  */
template <class TypeTag, class Assembler, class LinearSolver>
class RichardsPrivarSwitchNewtonSolverImplementation<TypeTag, Assembler, LinearSolver, true> : public PriVarSwitchNewtonSolver<Assembler, LinearSolver, GetPropType<TypeTag, Properties::PrimaryVariableSwitch>>
{
    using PrimaryVariableSwitch = GetPropType<TypeTag, Properties::PrimaryVariableSwitch>;
    using ParentType = PriVarSwitchNewtonSolver<Assembler, LinearSolver, PrimaryVariableSwitch>;
public:
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
