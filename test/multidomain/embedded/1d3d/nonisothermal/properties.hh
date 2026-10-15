// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_PROPERTIES_HH
#define BENCHMARKS_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/solidenergy/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/multidomain/traits.hh>

#include "wellborecouplingmanager.hh"

#include "problem_rock.hh"
#include "problem_well.hh"
#include "spatialparams_rock.hh"
#include "spatialparams_well.hh"

#ifndef COUPLINGMODE
#define COUPLINGMODE Embedded1d3dCouplingMode::Surface
#endif

#ifndef COORDINATETYPE
#define COORDINATETYPE TensorProductCoordinates
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Rock { using InheritsFrom = std::tuple<SolidEnergy, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>

struct Grid<TypeTag, TTag::Rock> { using type = Dune::YaspGrid<3, Dune::COORDINATETYPE<double, 3>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Rock> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Rock> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Rock> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Rock> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Rock> { using type = RockProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Rock>
{
    using type = RockSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};


// Create new type tags
namespace TTag {
struct Well { using InheritsFrom = std::tuple<OnePNI, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Well> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Well> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Well> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Well> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Well> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Well> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Well> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Well> { using type = WellProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Well>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Well> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Well>
{
    using type = WellSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};

template<class Traits>
using TheCouplingManager = WellboreCouplingManager<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Rock> { using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::Well>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Rock> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Rock> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Well> { using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::Rock, TypeTag>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Well> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Well> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
