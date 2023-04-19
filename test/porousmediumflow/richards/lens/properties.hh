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
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
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
struct RichardsLensBox { using InheritsFrom = std::tuple<RichardsLens, BoxModel>; };
struct RichardsLensCC { using InheritsFrom = std::tuple<RichardsLens, CCTpfaModel>; };
struct RichardsLensCCMpfa { using InheritsFrom = std::tuple<RichardsLens, CCMpfaModel>; };
} // end namespace TTag

#ifndef GRIDTYPE
// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = Dune::YaspGrid<2>; };
#else
// Use GRIDTYPE from CMakeLists.txt
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = GRIDTYPE; };
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsLens> { using type = RichardsLensProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsLens>
{
    using type = RichardsLensSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsLens>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} // end namespace Dumux::Properties

#endif
