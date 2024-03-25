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

#include "spatialparams.hh"
#include "problem.hh"

// Specify the properties for the lens problem
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsLens { using InheritsFrom = std::tuple<Richards>; };
struct RichardsLensBox { using InheritsFrom = std::tuple<Richards, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLensBox> { using type = Dune::YaspGrid<2>; };


// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsLensBox> { using type = RichardsLensProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsLensBox>
{
    using type = RichardsLensSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsLensBox>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} // end namespace Dumux::Properties

#endif
