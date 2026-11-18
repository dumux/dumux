// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The properties for the DFG benchmark of the laminar flow around a cylinder (Navier-Stokes)
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_FLOW_AROUND_CYLINDER_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_FLOW_AROUND_CYLINDER_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM DFGChannelTestMomentumDiamond
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS DFGChannelTestMassTpfa
#endif

#ifndef NEW_PROBLEM_INTERFACE
#define NEW_PROBLEM_INTERFACE 0
#endif

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"
#include "problem_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DFGChannelTest {};
struct DFGChannelTestMomentumDiamond { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct DFGChannelTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct DFGChannelTestMassTpfa { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
struct DFGChannelTestMassBox { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, BoxModel>; };
struct DFGChannelTestMassDiamond { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, FaceCenteredDiamondModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DFGChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DFGChannelTest>
{
    using type = Dune::UGGrid<2>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{
#if NEW_PROBLEM_INTERFACE
    using type = DFGChannelTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>;
#else
    using type = DFGChannelTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>;
#endif
};

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{
#if NEW_PROBLEM_INTERFACE
    using type = DFGChannelTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>;
#else
    using type = DFGChannelTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>;
#endif
};

#if NEW_PROBLEM_INTERFACE
//! The grid variables with the CVFE grid variables cache (momentum and mass)
template<class TypeTag>
using DFGChannelTestCVFEGridVariables = Dumux::Experimental::GridVariables<
    GetPropType<TypeTag, Properties::GridGeometry>,
    Dumux::Experimental::CVFE::CVFEGridVariablesCache<
        Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<
            GetPropType<TypeTag, Properties::Problem>,
            Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>,
            Dumux::CVFE::LocalBasisInterpolationPointData<GetPropType<TypeTag, Properties::GridGeometry>>
        >,
        getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>()
    >
>;

template<class TypeTag>
struct GridVariables<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = DFGChannelTestCVFEGridVariables<TypeTag>; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::TYPETAG_MASS>
{ using type = DFGChannelTestCVFEGridVariables<TypeTag>; };

//! The CVFE variables for the momentum problem
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TYPETAG_MOMENTUM>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMomentumCVFEVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMomentumCVFEVariables<Traits>;
};
#endif

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DFGChannelTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
