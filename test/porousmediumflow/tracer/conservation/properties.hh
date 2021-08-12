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
 * \ingroup TracerTests
 * \brief properties for the tracer conservation test
 */
#ifndef DUMUX_TEST_TRACER_CONSERVATION_PROPERTIES_HH
#define DUMUX_TEST_TRACER_CONSERVATION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/2p/model.hh>

#include "spatialparams_tracer.hh"
#include "problem_tracer.hh"
#include "spatialparams_2p.hh"
#include "problem_2p.hh"

namespace Dumux {

namespace Properties {

////////////////////////////////////////
// tracer properties
////////////////////////////////////////

//Create new type tags
namespace TTag {
struct TracerConservationTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerConservationTestTpfa { using InheritsFrom = std::tuple<TracerConservationTest, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerConservationTest> { using type = Dune::YaspGrid<1>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerConservationTest> { using type = TracerConservationTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerConservationTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerConservationTestSpatialParams<GridGeometry, Scalar>;
};

//! A simple fluid system with one tracer component
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                                TracerFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    //! If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! No component is the main component
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    //! The number of components
    static constexpr int numComponents = 1;

    //! Human readable component name (for vtk output)
    static std::string componentName(int compIdx)
    { return "Cs137"; }

    //! Molar mass in kg/mol
    static Scalar molarMass(int compIdx)
    { return 0.1; }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D = getParam<Scalar>("Problem.BinaryDiffusionCoefficient");
        return D;
    }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerConservationTest> { using type = TracerFluidSystem<TypeTag>; };

////////////////////////////////////////
// 2p flow properties
////////////////////////////////////////

// Create new type tags
namespace TTag {
struct TwoPFlow { using InheritsFrom = std::tuple<TwoP>; };
struct TwoPFlowTpfa { using InheritsFrom = std::tuple<TwoPFlow, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPFlow> { using type = Dune::YaspGrid<1>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPFlow> { using type = TwoPFlowTestProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPFlow>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPFlowTestSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties
} // end namespace Dumux
#endif
