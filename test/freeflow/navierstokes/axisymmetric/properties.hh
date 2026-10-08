// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Properties of the axisymmetric Stokes test (momentum only, rotation about the x-axis)
 *
 * TYPETAG selects the momentum scheme; NEW_PROBLEM_INTERFACE selects the general grid variables
 * and the integral problem interface, which the hybrid schemes require.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_AXISYMMETRIC_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_AXISYMMETRIC_PROPERTIES_HH

#ifndef TYPETAG
#define TYPETAG AxisymmetricStokesPQ1Bubble
#endif

#ifndef NEW_PROBLEM_INTERFACE
#define NEW_PROBLEM_INTERFACE 0
#endif

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct AxisymmetricStokes {};
struct AxisymmetricStokesPQ1Bubble { using InheritsFrom = std::tuple<AxisymmetricStokes, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct AxisymmetricStokesPQ1BubbleHybrid { using InheritsFrom = std::tuple<AxisymmetricStokes, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct AxisymmetricStokesPQ2Hybrid { using InheritsFrom = std::tuple<AxisymmetricStokes, NavierStokesMomentumCVFE, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::AxisymmetricStokes>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::AxisymmetricStokes>
{ using type = Dune::UGGrid<2>; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::AxisymmetricStokesPQ1Bubble>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    struct GGTraits : public PQ1BubbleDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<1>; };
public:
    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::AxisymmetricStokesPQ1BubbleHybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BaseTraits = HybridPQ1BubbleCVFEGridGeometryTraits<
        PQ1BubbleDefaultGridGeometryTraits<GridView, PQ1BubbleMapperTraits<GridView, /*numCubeBubbles*/2>>
    >;
    struct GGTraits : public BaseTraits
    { using Extrusion = RotationalExtrusion<1>; };
public:
    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::AxisymmetricStokesPQ2Hybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    struct GGTraits : public PQ2DefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<1>; };
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::AxisymmetricStokes>
{
#if NEW_PROBLEM_INTERFACE
    using type = AxisymmetricStokesTestProblem<TypeTag, CVFENavierStokesMomentumProblem<TypeTag>>;
#else
    using type = AxisymmetricStokesTestProblem<TypeTag, NavierStokesMomentumProblem<TypeTag>>;
#endif
};

#if NEW_PROBLEM_INTERFACE
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::AxisymmetricStokes>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = NavierStokesMomentumCVFEVariables<NavierStokesMomentumCVFEVolumeVariablesTraits<PV, FSY, FST, MT>>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::AxisymmetricStokesPQ1Bubble>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};
#endif

} // end namespace Dumux::Properties

#endif
