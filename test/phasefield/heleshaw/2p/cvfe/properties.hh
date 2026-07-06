// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawTests
 * \brief The properties of a test problem for the hybrid CVFE Hele-Shaw
 *        two-phase model, discretized with quadratic (PQ2 hybrid CVFE)
 *        Lagrange elements for pressure, phase field, and chemical potential.
 */
#ifndef DUMUX_TEST_PHASEFIELD_HELESHAW_2P_CVFE_PROPERTIES_HH
#define DUMUX_TEST_PHASEFIELD_HELESHAW_2P_CVFE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/phasefield/heleshaw/2p/cvfe/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct TestHeleShawTwoPCVFE { using InheritsFrom = std::tuple<HeleShawTwoPCVFE, PQ2HybridModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestHeleShawTwoPCVFE> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestHeleShawTwoPCVFE> { using type = Dumux::HeleShawTwoPCVFETestProblem<TypeTag>; };

// PQ2FVGridGeometry reads this property, so the geometry cache must be enabled.
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TestHeleShawTwoPCVFE> : public std::true_type {};

// The grid variables: PQ2 always uses the hybrid CVFE local dof/quadrature machinery,
// bundling variables and interpolation-point data caches (gradients, shape values)
// in one grid variables cache.
template<class TypeTag>
struct GridVariables<TypeTag, TTag::TestHeleShawTwoPCVFE>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GridGeometry>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GridVariablesCache = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, /*enableCache*/true>;
public:
    using type = Dumux::Experimental::GridVariables<GridGeometry, GridVariablesCache>;
};

} // end namespace Dumux::Properties

#endif
