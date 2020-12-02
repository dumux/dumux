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
#ifndef DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_IMPES_ADAPTIVE_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

#include "test_impesadaptiveproblem.hh"
#include "test_impesadaptivespatialparams.hh"

#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct TestIMPESAdaptive { using InheritsFrom = std::tuple<TestIMPESAdaptiveSpatialParams, IMPESTwoPAdaptive, FVTransportTwoP, FVPressureTwoPAdaptive>; };
struct TestIMPESAdaptiveRestart { using InheritsFrom = std::tuple<TestIMPESAdaptive>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TestIMPESAdaptive> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestIMPESAdaptive> { using type = TestIMPESAdaptiveProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestIMPESAdaptive>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};
} // namespace Properties

#endif
