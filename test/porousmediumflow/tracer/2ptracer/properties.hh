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
 * \brief The properties of the 2p problem with multiple tracer bands in a porous groundwater reservoir with a lens
 */
#ifndef DUMUX_2P_TRACER_TEST_PROPERTIES_HH
#define DUMUX_2P_TRACER_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <test/porousmediumflow/2p/incompressible/properties.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

//Create new type tags
namespace TTag {
struct TwoPTracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TwoPTracerTestTpfa { using InheritsFrom = std::tuple<TwoPTracerTest, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPTracerTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPTracerTest> { using type = TwoPTracerTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPTracerTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTracerTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPTracerTest> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::TwoPTracerTestTpfa> { static constexpr bool value = false; };

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

    //! The number of components
    static constexpr int numComponents = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

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
struct FluidSystem<TypeTag, TTag::TwoPTracerTest> { using type = TracerFluidSystem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
