// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the single-phase flow
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
#ifndef DUMUX_GEOMECHANICS_ELONEP_PROPERTIES_HH
#define DUMUX_GEOMECHANICS_ELONEP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

#include "spatialparams_1p.hh"
#include "spatialparams_poroelastic.hh"

#include "problem_1p.hh"
#include "problem_poroelastic.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePSub { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSub>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSub> { using type = Dune::YaspGrid<2>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSub> { using type = OnePSubProblem<TypeTag> ; };
// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSub>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using type = OnePSpatialParams<GridGeometry, Scalar, CouplingManager>;
};

// Create new type tags
namespace TTag {
struct PoroElasticSub { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = Dune::YaspGrid<2>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticSub> { using type = Dumux::PoroElasticSubProblem<TypeTag>; };
// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticSub>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};
// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticSub>
{
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry>,
                                           CouplingManager>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<TTag::OnePSub, TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

} // end namespace Dumux::Properties

#endif
