// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief The properties for the  one-phase two-component pore network model.
 */
#ifndef DUMUX_PNM1P2C_PROPERTIES_HH
#define DUMUX_PNM1P2C_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porenetwork/1pnc/model.hh>
#include <dumux/porenetwork/1p/spatialparams.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include "problem.hh"
#include "spatialparams.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNC>; };
#else
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNCNI>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePTwoCProblem> { using type = PNMOnePTwoCProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePTwoCProblem>
{
    using Policy = FluidSystems::H2ON2DefaultPolicy</*simple*/true>;
    using H2ON2 = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, Policy>;
    using type = FluidSystems::OnePAdapter<H2ON2>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOnePTwoCProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::CompositionalSpatialParams<GridGeometry, Scalar>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePTwoCProblem> { using type = Dune::FoamGrid<1, 3>; };

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePTwoCProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityAzzamDullien<Scalar>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

} //end namespace Dumux::Properties

#endif
