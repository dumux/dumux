// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerTests
 * \brief The properties for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
#ifndef DUMUX_TRACER_TEST_PROPERTIES_HH
#define DUMUX_TRACER_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams.hh"

#ifndef USEMOLES // default to true if not set through CMake
#define USEMOLES true
#endif

#ifndef ENABLEDISPERSION // default to false if not set through CMake
#define ENABLEDISPERSION false
#endif

#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestTpfa { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
struct TracerTestMpfa { using InheritsFrom = std::tuple<TracerTest, CCMpfaModel>; };
struct TracerTestBox { using InheritsFrom = std::tuple<TracerTest, BoxModel>; };
} // end namespace TTag

// enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTest<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = USEMOLES; };

template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::TracerTest> { static constexpr bool value = ENABLEDISPERSION; };

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
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! The number of components
    static constexpr int numComponents = 2;
    static constexpr int numPhases = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }

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
        static const Scalar D = getParam<Scalar>("Problem.D");
        static const Scalar D2 = getParam<Scalar>("Problem.D2");
        if (compIdx == 0)
            return D;
        else
            return D2;
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
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = TracerFluidSystem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
