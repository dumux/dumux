// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Properties for the Stokes Darcy problem
 */

#ifndef DUMUX_DARCYSTOKES_PROPERTIES_HH
#define DUMUX_DARCYSTOKES_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#include <dumux/porousmediumflow/2p/model.hh>

#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>
#include <dumux/multidomain/staggeredtraits.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyTwoP { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyTwoP> { using type = Dumux::DarcySubProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyTwoP> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyTwoP> { static constexpr bool value = false; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyTwoP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConservationSpatialParams<GridGeometry, Scalar>;
};

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::DarcyTwoP>
{ static constexpr auto value = TwoPFormulation::p1s0; };

// Create new type tags
namespace TTag {
struct StokesOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOneP>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyTwoP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyTwoP>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOneP, Properties::TTag::StokesOneP, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingFluidSystem
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OType = Dumux::Components::SimpleH2O<Scalar>;
    using H2OPhase = Dumux::FluidSystems::OnePLiquid<Scalar, H2OType>;
    using AirType = Dumux::Components::TabulatedComponent<Components::Air<Scalar>, false >;
    using AirPhase = Dumux::FluidSystems::OnePGas<Scalar, AirType>;
    using type = FluidSystems::TwoPImmiscible<Scalar, H2OPhase, AirPhase> ;
};

// the fluid system for the free-flow model
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP>
{
    static constexpr auto phaseIdx = 1; // simulate the air phase
    using type = FluidSystems::OnePAdapter<typename CouplingFluidSystem<TypeTag>::type, phaseIdx>;
};

// the fluid system for the Darcy model
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyTwoP>
{
    using type = typename CouplingFluidSystem<TypeTag>::type;
};

} // end namespace Dumux::Properties

#endif
