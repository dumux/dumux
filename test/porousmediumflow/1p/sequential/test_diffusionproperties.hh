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
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
#ifndef DUMUX_TEST_2P_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_2P_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/omethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>

#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

#include "test_diffusionproblem.hh"
#include "test_diffusionspatialparams.hh"

namespace Dumux::Properties
{
//// set the types for the 2PFA FV method
// Create new type tags
namespace TTag {

struct FVVelocity2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, FVPressureTwoP>; };

// set the types for the MPFA-O FV method
struct FVMPFAOVelocity2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, FvMpfaO2dPressureTwoP>; };

// set the types for the mimetic FD method
struct MimeticPressure2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, MimeticPressureTwoP>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::FVVelocity2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FVVelocity2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FVVelocity2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

//template<class TypeTag>
// struct LinearSolver<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = ILUnBiCGSTABBackend; };
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = SSORBiCGSTABBackend; };
template<class TypeTag>
struct Problem<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};


template<class TypeTag>
struct Problem<TypeTag, TTag::MimeticPressure2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MimeticPressure2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MimeticPressure2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

} // end namespace Properties

#endif
