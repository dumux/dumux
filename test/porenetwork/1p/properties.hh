// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief The properties for the one-phase pore network model.
 */
#ifndef DUMUX_PNM1P_PROPERTIES_HH
#define DUMUX_PNM1P_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>

// Pore network model
#include <dumux/porenetwork/1p/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// the problem
#include "problem.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

namespace TTag {
#if ISOTHERMAL
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOneP>; };
#else
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOnePNI>; };
#endif
}

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePProblem>
{ using type = PNMOnePProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePProblem>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> > ;
};

// the grid
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePProblem>
{ using type = Dune::FoamGrid<1, 3>; };

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type =  PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMOnePProblem> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

} //end namespace Dumux::Properties

#endif
