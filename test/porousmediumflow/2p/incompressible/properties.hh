// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p test.
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef ENABLEINTERFACESOLVER
#define ENABLEINTERFACESOLVER 0
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPIncompressible { using InheritsFrom = std::tuple<TwoP>; };
struct TwoPIncompressibleTpfa { using InheritsFrom = std::tuple<TwoPIncompressible, CCTpfaModel>; };
struct TwoPIncompressibleMpfa { using InheritsFrom = std::tuple<TwoPIncompressible, CCMpfaModel>; };
struct TwoPIncompressibleBox { using InheritsFrom = std::tuple<TwoPIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressible> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressible> { using type = TwoPTestProblem<TypeTag>; };

// the local residual containing the analytic derivative methods
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPIncompressible> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPIncompressible>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPTestSpatialParams<GridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = false; };

// Maybe enable the box-interface solver
template<class TypeTag>
struct EnableBoxInterfaceSolver<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = ENABLEINTERFACESOLVER; };

} // end namespace Dumux::Properties

#endif
