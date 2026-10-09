// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Properties of single-phase flow through the pore network of a sphere packing
 */
#ifndef DUMUX_PNM_EXTRACTION_PERMEABILITY_PROPERTIES_HH
#define DUMUX_PNM_EXTRACTION_PERMEABILITY_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porenetwork/1p/model.hh>
#include <dumux/porenetwork/extraction/spatialparams.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct SpherePackingPermeability { using InheritsFrom = std::tuple<PNMOneP>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::SpherePackingPermeability>
{ using type = SpherePackingPermeabilityProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SpherePackingPermeability>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::SpherePackingPermeability>
{ using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SpherePackingPermeability>
{
    using type = PoreNetwork::SpherePackingOnePSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                            GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::SpherePackingPermeability>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = PoreNetwork::CreepingFlow<Scalar, PoreNetwork::TransmissibilityChareyre<Scalar>>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SpherePackingPermeability>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
