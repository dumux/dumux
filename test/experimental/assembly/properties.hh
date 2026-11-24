// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The properties for the analytic-Jacobian assembly test
 */
#ifndef DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_PROPERTIES_HH
#define DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/concepts/ipdata_.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/volumevariables.hh>

#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>

#include "localresidual.hh"
#include "problem.hh"

namespace Dumux {

struct NonlinearDiffusionIndices
{ static constexpr int equationIdx(int i) { return i; } };

struct NonlinearDiffusionModelTraits
{
    using Indices = NonlinearDiffusionIndices;
    static constexpr int numEq() { return 1; }
};

/*!
 * \brief The primary variables at a local dof, and nothing else
 *
 * Deliberately *not* `BasicVolumeVariables`: that one updates itself from a sub-control volume,
 * and `Detail::CVFE::VariablesAdapter` then dispatches through `fvGeometry.scv(localDofIndex)`,
 * which a hybrid scheme's Galerkin-tested dofs do not have. Offering only the interpolation-point
 * interface is what makes a model usable on every CVFE scheme.
 */
template<class Traits>
class NonlinearDiffusionVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    //! `Concept::LocalDofIpData` is load-bearing, not decoration: an unconstrained template
    //! also matches `update(elemSol, problem, element, scv)`, and `VariablesAdapter` detects
    //! that signature by whether it compiles. It would then dispatch through
    //! `fvGeometry.scv(localDofIndex)`, which is out of range for a Galerkin-tested local dof.
    template<class ElementSolution, class Problem, class ElementDiscretization,
             Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const ElementDiscretization& elemDisc,
                const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    Scalar priVar(const int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    Scalar extrusionFactor() const { return 1.0; }

private:
    PrimaryVariables priVars_;
};

template<class PV, class MT>
struct NonlinearDiffusionVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct NonlinearDiffusion { using InheritsFrom = std::tuple<ModelProperties>; };
struct NonlinearDiffusionPQ1Bubble { using InheritsFrom = std::tuple<NonlinearDiffusion, PQ1BubbleModel>; };
struct NonlinearDiffusionPQ1BubbleHybrid { using InheritsFrom = std::tuple<NonlinearDiffusion, PQ1BubbleHybridModel>; };
struct NonlinearDiffusionPQ2Hybrid { using InheritsFrom = std::tuple<NonlinearDiffusion, PQ2HybridModel>; };
struct NonlinearDiffusionPQ3Hybrid { using InheritsFrom = std::tuple<NonlinearDiffusion, PQ3HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::NonlinearDiffusion>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::NonlinearDiffusion>
{ using type = NonlinearDiffusionTestProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NonlinearDiffusion>
{ using type = NonlinearDiffusionLocalResidual<TypeTag>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NonlinearDiffusion>
{ using type = NonlinearDiffusionModelTraits; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NonlinearDiffusion>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = NonlinearDiffusionVolumeVariables<NonlinearDiffusionVolumeVariablesTraits<PV, MT>>;
};

//! selecting the general grid variables is what puts this model on the experimental stack
template<class TypeTag>
struct GridVariables<TypeTag, TTag::NonlinearDiffusion>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableHybrid = []{
        if constexpr (requires { GG::enableHybridCVFE; })
            return GG::enableHybridCVFE;
        else
            return false;
    }();
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = std::conditional_t<enableHybrid,
        Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>,
        Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>>;
    using GVC = std::conditional_t<enableHybrid,
        Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>,
        Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

//! the general grid variables cache only implements the caching-enabled specialisation
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::NonlinearDiffusion>
{ static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
