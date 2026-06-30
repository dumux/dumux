// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_PROPERTIES_HH
#define DUMUX_COOKS_MEMBRANE_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/elastic/model.hh>

// New variables concept
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct CooksMembrane
{
    using InheritsFrom = std::tuple<Elastic, DISCRETIZATION_MODEL>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembrane>
{ using type = CooksMembraneProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembrane>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneSpatialParams<Scalar, GridGeometry>;
};

//! Use Experimental::GridVariables for all three discretizations so that
//! ElementVariables covers all local DOFs via the appropriate GVC.
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembrane>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Prob = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using DiscMethod = typename GG::DiscretizationMethod;

    // PQ2 uses HybridCVFEGridVariablesCache; Box and PQ1Bubble use the standard CVFEGridVariablesCache_.
    static constexpr bool isPQ2 = std::is_same_v<DiscMethod, DiscretizationMethods::PQ2>;

    using HybridTraits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Prob, Variables, IPDataCache>;
    using StandardTraits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Prob, Variables, IPDataCache>;

    using GVC = std::conditional_t<isPQ2,
        Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<HybridTraits, enableCache>,
        Dumux::Experimental::CVFE::CVFEGridVariablesCache<StandardTraits, enableCache>>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

} // end namespace Dumux::Properties

#endif
