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
 * \ingroup SequentialTwoPTests
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_IMPES_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_IMPES_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

//following includes are only needed if a global pressure formulation is chosen!
//Then only a total velocity can be reconstructed for the transport step
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/capillarydiffusion.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/gravitypart.hh>

#include "test_impesproblem.hh"
#include "test_impesspatialparams.hh"

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct IMPESTest { using InheritsFrom = std::tuple<TestIMPESSpatialParams, FVTransportTwoP, IMPESTwoP, FVPressureTwoP>; };

// set up an additional problem where the AMG backend is used
struct IMPESTestWithAMG { using InheritsFrom = std::tuple<IMPESTest>; };

} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::IMPESTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::IMPESTest> { using type = IMPESTestProblem<TypeTag>; };

////////////////////////////////////////////////////////////////////////
//Switch to a p_n-S_w formulation
//
// template<class TypeTag>
// struct Formulation<TypeTag, TTag::IMPESTest> { static constexpr int value = SequentialTwoPCommonIndices::pnsn; };
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Switch to a p_global-S_w formulation
//
// template<class TypeTag>
// struct Formulation<TypeTag, TTag::IMPESTest> { static constexpr int value = SequentialTwoPCommonIndices::pGlobalSw; };
//
//Define the capillary pressure term in the transport equation -> only needed in case of a p_global-S_w formulation!
template<class TypeTag>
struct CapillaryFlux<TypeTag, TTag::IMPESTest> { using type = CapillaryDiffusion<TypeTag>; };
//
//Define the gravity term in the transport equation -> only needed in case of a p_global-S_w formulation!
template<class TypeTag>
struct GravityFlux<TypeTag, TTag::IMPESTest> { using type = GravityPart<TypeTag>; };
//
////////////////////////////////////////////////////////////////////////

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IMPESTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::IMPESTest> { using type = EvalCflFluxCoats<TypeTag>; };

// use the AMG backend for the corresponding test
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::IMPESTestWithAMG>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;

};
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::IMPESTestWithAMG> { using type = Dune::YaspGrid<2>; };

} // end namespace Properties


#endif
