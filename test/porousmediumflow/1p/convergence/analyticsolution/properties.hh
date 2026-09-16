// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the convergence test with analytic solution
 */
#ifndef DUMUX_CONVERGENCE_TEST_ONEP_PROPERTIES_HH
#define DUMUX_CONVERGENCE_TEST_ONEP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq3.hh>

#include <dumux/flux/cvfe/darcyslaw_.hh>
#include <dumux/porousmediumflow/1p/variables.hh>
#include <dumux/porousmediumflow/immiscible/localresidual_.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"
#include "problem.hh"
#include "problem_newinterface.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePConvergence { using InheritsFrom = std::tuple<OneP>; };
struct OnePConvergenceTpfa { using InheritsFrom = std::tuple<OnePConvergence, CCTpfaModel>; };
struct OnePConvergenceMpfa { using InheritsFrom = std::tuple<OnePConvergence, CCMpfaModel>; };
struct OnePConvergenceBox { using InheritsFrom = std::tuple<OnePConvergence, BoxModel>; };

//! Everything the schemes assembled with the new interfaces have in common
struct OnePConvergenceNewInterface { using InheritsFrom = std::tuple<OnePConvergence>; };
struct OnePConvergencePQ1Bubble { using InheritsFrom = std::tuple<OnePConvergenceNewInterface, PQ1BubbleHybridModel>; };
struct OnePConvergencePQ3 { using InheritsFrom = std::tuple<OnePConvergenceNewInterface, PQ3HybridModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePConvergence> { using type = Dumux::ConvergenceProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePConvergence>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePConvergence> { using type = GRIDTYPE; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePConvergence>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };

//! The problem stated in terms of the new interfaces
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePConvergenceNewInterface>
{ using type = Dumux::ConvergenceProblemNewInterface<TypeTag>; };

/*!
 * \brief The grid variables holding the quantities per local dof
 * \note The variables are asked for at an interpolation point, which is what the degrees
 *       of freedom owning no sub-control volume require.
 */
template<class TypeTag>
struct GridVariables<TypeTag, TTag::OnePConvergenceNewInterface>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using VariablesTraits = OnePVolumeVariablesTraits<
        GetPropType<TypeTag, Properties::PrimaryVariables>,
        GetPropType<TypeTag, Properties::FluidSystem>,
        GetPropType<TypeTag, Properties::FluidState>,
        GetPropType<TypeTag, Properties::SolidSystem>,
        GetPropType<TypeTag, Properties::SolidState>,
        typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType,
        GetPropType<TypeTag, Properties::ModelTraits>
    >;
    using Variables = Dumux::Experimental::OnePVariables<VariablesTraits>;

    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;

public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

//! Evaluate Darcy's law at the interpolation points
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePConvergenceNewInterface>
{
    using type = Experimental::CVFEDarcysLawAtIp<GetPropType<TypeTag, Properties::Scalar>,
                                                 GetPropType<TypeTag, Properties::GridGeometry>>;
};

//! Assemble with the local residual that leaves the quadrature to the scheme
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePConvergenceNewInterface>
{ using type = Experimental::ImmiscibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
