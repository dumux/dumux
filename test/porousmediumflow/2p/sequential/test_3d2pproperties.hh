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
 * \brief test problem for sequential 2p models in 3d
 */

#ifndef DUMUX_TEST_3D2P_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_3D2P_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressurepropertiesadaptive.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

#include "test_3d2pproblem.hh"
#include "test_3d2pspatialparams.hh"

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/gridadaptionindicatorlocal.hh>


namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct ThreeDTwoPTest { using InheritsFrom = std::tuple<Test3d2pSpatialParams>; };
struct FVTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, FVPressureTwoP>; };
struct FVAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FVPressureTwoPAdaptive>; };
struct MPFALTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaL3dPressureTwoP>; };
struct MPFALAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FvMpfaL3dPressureTwoPAdaptive>; };
struct MimeticTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, MimeticPressureTwoP>; };
struct MimeticAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, MimeticPressureTwoPAdaptive>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreeDTwoPTest> { using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreeDTwoPTest> { using type = Test3D2PProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreeDTwoPTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

#if PROBLEM == 1
template<class TypeTag>
struct Formulation<TypeTag, TTag::ThreeDTwoPTest> { static constexpr int value = SequentialTwoPCommonIndices::pnSw; };
template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::ThreeDTwoPTest> { using type = EvalCflFluxCoats<TypeTag>; };
#endif

template<class TypeTag>
struct AdaptionIndicator<TypeTag, TTag::ThreeDTwoPTest> { using type = GridAdaptionIndicator2PLocal<TypeTag>; };
} // end namespace Properties

#endif
