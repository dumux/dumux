// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCTests
 * \brief The properties of the problem showcasing the capabilities of the kinetic model.
 */
#ifndef DUMUX_EVAPORATION_ATMOSPHERE_PROPERTIES_HH
#define DUMUX_EVAPORATION_ATMOSPHERE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/h2on2kinetic.hh>
#include <dumux/porousmediumflow/mpnc/pressureformulation.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct EvaporationAtmosphere { using InheritsFrom = std::tuple<MPNCNonequil>; };
struct EvaporationAtmosphereBox { using InheritsFrom = std::tuple<EvaporationAtmosphere, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::EvaporationAtmosphere>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::EvaporationAtmosphere>
{ using type = EvaporationAtmosphereProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::EvaporationAtmosphere>
{
    using type = FluidSystems::H2ON2Kinetic<GetPropType<TypeTag, Properties::Scalar>,
                                            FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

//! Set the default pressure formulation: either pw first or pn first
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::EvaporationAtmosphere>
{
public:
    static const MpNcPressureFormulation value = MpNcPressureFormulation::leastWettingFirst;
};

// Set the fluid system
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::EvaporationAtmosphere>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::EvaporationAtmosphere>
{
    using GridGeometry = GetPropType<TypeTag, GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = EvaporationAtmosphereSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
