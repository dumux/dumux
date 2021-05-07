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
 * \brief The properties for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
#ifndef DUMUX_TRACER_TEST_PROPERTIES_HH
#define DUMUX_TRACER_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/freeflow/advectiondiffusion/model.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams.hh"
#include "problem_flow.hh"
#include "problem_transport.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct NavierStokesTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
struct AdvectionDiffusionTest { using InheritsFrom = std::tuple<AdvectionDiffusion>; };
struct AdvectionDiffusionTpfa { using InheritsFrom = std::tuple<AdvectionDiffusionTest, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::NavierStokesTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::NavierStokesTest>
{
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim> >;
    using HostGrid = TensorGrid;
    using type = Dune::SubGrid<dim, HostGrid>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::NavierStokesTest> { using type = Dumux::FlowProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::NavierStokesTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::NavierStokesTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::NavierStokesTest> { static constexpr bool value = true; };

/////////////////////////////////////
// The Transport Problem Properties
/////////////////////////////////////

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::AdvectionDiffusionTest>
{
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim> >;
    using HostGrid = TensorGrid;
    using type = Dune::SubGrid<dim, HostGrid>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::AdvectionDiffusionTest> { using type = AdvectionDiffusionTest<TypeTag>; };

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::AdvectionDiffusionTest> { using type = PorousMediumFluxVariablesCache<TypeTag>; };
template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::AdvectionDiffusionTest> { using type = PorousMediumFluxVariablesCacheFiller<TypeTag>; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::AdvectionDiffusionTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = AdvectionDiffusionTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::AdvectionDiffusionTest> { static constexpr bool value = true; };

//! A simple fluid system with one tracer component
template<class TypeTag>
class TransportFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                    TransportFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! The number of components
    static constexpr int numComponents = 1;
    static constexpr int numPhases = 1;

//     //! Molar mass in kg/mol of the component with index compIdx
//     static Scalar molarMass(unsigned int compIdx)
//     { return 18; }
//
//     static Scalar density(unsigned int compIdx)
//     { return 1000.0; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "TransportedComponent_" + std::to_string(compIdx); }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D1 = getParam<Scalar>("Problem.D1");
        if (compIdx == 0)
            return D1;
    }

    /*!
     * \copydoc Dumux::FluidSystems::Base::isCompressible
     */
    static constexpr bool isCompressible(int phaseIdx)
    { return false; }

     /*!
     * \copydoc  Dumux::FluidSystems::Base::viscosityIsConstant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return true; }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::AdvectionDiffusionTest> { using type = TransportFluidSystem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
