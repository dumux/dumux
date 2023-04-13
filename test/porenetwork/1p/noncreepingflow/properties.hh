// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief The properties for the one-phase pore-network model with non-creeping flow.
 */
#ifndef DUMUX_PNM1P_NONCREEPING_PROPERTIES_HH
#define DUMUX_PNM1P_NONCREEPING_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/1p/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"
#include "../properties.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

namespace TTag {
struct PNMOnePNonCreepingProblem { using InheritsFrom = std::tuple<PNMOnePProblem>; };
}

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOnePNonCreepingProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::NonCreepingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNonCreepingProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw =  PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::NonCreepingFlow<Scalar, TransmissibilityLaw>;
};

} //end namespace Dumux::Properties

#endif
