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
 * \brief test problem for sequential 2p models
 */

#ifndef DUMUX_TEST_MPFA2P_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_MPFA2P_PROBLEM_PROPERTIES_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/trichloroethene.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/omethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/2dpressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/gridadaptionindicatorlocal.hh>

#include "test_mpfa2pproblem.hh"
#include "test_mpfa2pspatialparams.hh"

namespace Dumux::Properties
{

// Create new type tags
namespace TTag {
struct MPFATwoPTest { using InheritsFrom = std::tuple<Test2PSpatialParams>; };
struct FVTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FVPressureTwoP>; };
struct FVAdaptiveTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FVPressureTwoPAdaptive>; };
struct MPFAOTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaO2dPressureTwoP>; };
struct MPFALTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaL2dPressureTwoP>; };
struct MPFALAdaptiveTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FvMpfaL2dPressureTwoPAdaptive>; };

} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::MPFATwoPTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::MPFATwoPTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MPFATwoPTest> { using type = MPFATwoPTestProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MPFATwoPTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#if PROBLEM == 2
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
#else
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#endif
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

#if PROBLEM == 1
template<class TypeTag>
struct Formulation<TypeTag, TTag::MPFATwoPTest> { static constexpr int value = SequentialTwoPCommonIndices::pnsw; };
#endif

template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::MPFATwoPTest> { using type = EvalCflFluxCoats<TypeTag>; };
template<class TypeTag>
struct AdaptionIndicator<TypeTag, TTag::MPFATwoPTest> { using type = GridAdaptionIndicator2PLocal<TypeTag>; };
} // end namespace Properties
#endif
