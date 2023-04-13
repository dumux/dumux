// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the two-phase flow
 *        sub-problem in the coupled poro-mechanical elp problem.
 */

#ifndef DUMUX_GEOMECHANICS_ELTWOP_PROPERTIES_HH
#define DUMUX_GEOMECHANICS_ELTWOP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/components/co2.hh>
#include <dumux/material/components/defaultco2table.hh>
#include <dumux/material/fluidsystems/brineco2.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/geomechanics/poroelastic/couplingmanager.hh>

#include "spatialparams_2p.hh"
#include "spatialparams_poroelastic.hh"

#include "problem_2p.hh"
#include "problem_poroelastic.hh"


namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPSub { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
} // end namespace TTag

// Set the fluid system for TwoPSubProblem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineCO2<Scalar, Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPSub> { using type = Dune::YaspGrid<3>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPSub> { using type = TwoPSubProblem<TypeTag> ; };
// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPSub>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using type = TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>;
};

// Create new type tags
namespace TTag {
struct PoroElasticSub { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticSub> { using type = Dune::YaspGrid<3>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticSub> { using type = Dumux::PoroElasticSubProblem<TypeTag>; };

// Set the fluid system for TwoPSubProblem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineCO2<Scalar, Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>>;
};

// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticSub>
{
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry>,
                                           CouplingManager, FluidSystem>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TwoPSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<Properties::TTag::TwoPSub, Properties::TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticSub>
{
private:
    // define traits etc. as below in main
    using Traits = MultiDomainTraits<Properties::TTag::TwoPSub, Properties::TTag::PoroElasticSub>;
public:
    using type = PoroMechanicsCouplingManager< Traits >;
};

} // end namespace Dumux::Properties

#endif
