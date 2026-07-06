// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnTests
 * \brief The properties of a test problem for the Allen-Cahn model
 *        discretized with quadratic (PQ2 hybrid CVFE) Lagrange elements.
 */
#ifndef DUMUX_TEST_PHASEFIELD_ALLENCAHN_PROPERTIES_HH
#define DUMUX_TEST_PHASEFIELD_ALLENCAHN_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/phasefield/allencahn/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct TestAllenCahnPQ2 { using InheritsFrom = std::tuple<AllenCahn, PQ2HybridModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestAllenCahnPQ2> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestAllenCahnPQ2> { using type = Dumux::AllenCahnTestProblem<TypeTag>; };

// PQ2FVGridGeometry reads this property, so the geometry cache must be enabled.
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TestAllenCahnPQ2> : public std::true_type {};

// The double-well reaction term f'(phi) is cubic in phi (phi itself
// degree-2 on PQ2), so as a function of position it is degree 6; combined
// with a degree-2 test function the weak-form reaction integral needs an
// exact degree-8 rule -- well above PQ2's default element-quadrature order
// (4, tuned for the linear flux terms). The grid variables cache is sized
// to the grid geometry's OWN quadrature scheme, so this must be raised here
// (not via an ad-hoc order override at the call site, which would read past
// the cache).
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TestAllenCahnPQ2>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using QuadratureTraits = PQ2QuadratureTraits<GridView,
        QuadratureRules::MidpointQuadrature, QuadratureRules::DuneQuadrature<2>,
        QuadratureRules::DuneQuadrature<8>, QuadratureRules::DuneQuadrature<4>, QuadratureRules::DuneQuadrature<4>>;
    using GGTraits = PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadratureTraits>;
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// The grid variables: PQ2 always uses the hybrid CVFE local dof/quadrature machinery,
// bundling variables and interpolation-point data caches (gradients, shape values)
// in one grid variables cache (see dumux/freeflow/navierstokes/momentum/cvfe/ for the
// same pattern used by the PQ2Hybrid Navier-Stokes momentum tests).
template<class TypeTag>
struct GridVariables<TypeTag, TTag::TestAllenCahnPQ2>
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
