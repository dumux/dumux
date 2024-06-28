// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief The properties of the water infiltration problem with a
 *        low-permeability lens embedded into a high-permeability domain which
 *        uses the Richards box model.
 */
#ifndef DUMUX_RICHARDS_LENSPROPERTIES_HH
#define DUMUX_RICHARDS_LENSPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/droplet/volumevariables.hh>
#include <dumux/porousmediumflow/droplet/iofields.hh>

#include "spatialparams.hh"
#include "problem.hh"

// Specify the properties for the lens problem
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsBoxDroplet { using InheritsFrom = std::tuple<Richards, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsBoxDroplet> { using type = Dune::SubGrid<3, Dune::YaspGrid<3>>; };


// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsBoxDroplet> { using type = RichardsProblemDroplet<TypeTag>; };


//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RichardsBoxDroplet>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = RichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = RichardsVolumeVariablesDroplet<Traits>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsBoxDroplet>
{
    using type = RichardsTestSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsBoxDroplet>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

//! Set the vtk output fields specific to the twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::RichardsBoxDroplet> { using type = TwoPIOFieldsDroplet; };
} // end namespace Dumux::Properties

#endif
