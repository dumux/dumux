// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The properties of the problem where air is injected under a low permeable layer in a depth of 2700m.
 */

#ifndef DUMUX_INJECTION_PROPERTIES_HH
#define DUMUX_INJECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/staticinteractionvolume.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Injection { using InheritsFrom = std::tuple<TwoPTwoC>; };
struct InjectionBox { using InheritsFrom = std::tuple<Injection, BoxModel>; };
struct InjectionCCTpfa { using InheritsFrom = std::tuple<Injection, CCTpfaModel>; };
struct InjectionCCMpfa { using InheritsFrom = std::tuple<Injection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection> { using type = InjectionProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection>
{
    using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Injection> { static constexpr bool value = true; };

// Enable caching or not (reference solutions created without caching)
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };

// use the static interaction volume around interior vertices in the mpfa test
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::InjectionCCMpfa>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // structured two-d grid
    static constexpr int numIvScvs = 4;
    static constexpr int numIvScvfs = 4;

    // use the default traits
    using Traits = CCMpfaODefaultStaticInteractionVolumeTraits< NodalIndexSet, Scalar, numIvScvs, numIvScvfs >;
public:
    using type = CCMpfaOStaticInteractionVolume< Traits >;
};

} // end namespace Dumux::Properties

#endif
