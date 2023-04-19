// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Properties of non-isothermal steam injection test problem for the 2p1cni model.
 */

#ifndef DUMUX_STEAM_INJECTIONPROBLEM_PROPERTIES_HH
#define DUMUX_STEAM_INJECTIONPROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2p1c/model.hh>

#include <dumux/material/fluidsystems/2p1c.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct InjectionProblem { using InheritsFrom = std::tuple<TwoPOneCNI>; };
struct TwoPOneCNIBox { using InheritsFrom = std::tuple<InjectionProblem, BoxModel>; };
struct TwoPOneCNICCTpfa { using InheritsFrom = std::tuple<InjectionProblem, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::InjectionProblem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InjectionProblem> { using type = InjectionProblem<TypeTag>; };


// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InjectionProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2OType = Dumux::Components::TabulatedComponent<Dumux::Components::H2O<Scalar> >;
public:
    using type = Dumux::FluidSystems::TwoPOneC<Scalar, H2OType >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InjectionProblem>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionProblemSpatialParams<GridGeometry, Scalar>;
};

//Define whether spurious cold-water flow into the steam is blocked
template<class TypeTag>
struct UseBlockingOfSpuriousFlow<TypeTag, TTag::InjectionProblem> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
