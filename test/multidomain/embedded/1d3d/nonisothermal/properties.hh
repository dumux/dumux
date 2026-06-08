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
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_soil.hh"
#include "problem_voids.hh"
#include "spatialparams_soil.hh"
#include "spatialparams_voids.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Soil { using InheritsFrom = std::tuple<OnePNI, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>

struct Grid<TypeTag, TTag::Soil> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Soil> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Soil> { using type = SoilProblem<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Soil> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Soil>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Water = Components::SimpleH2O<Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Water >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Soil>
{
    using type = SoilSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};


// Create new type tags
namespace TTag {
struct Voids { using InheritsFrom = std::tuple<OnePNI, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Voids> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Voids> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Voids> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Voids> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Voids> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Voids> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Voids> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Voids> { using type = VoidsProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Voids>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Water = Components::SimpleH2O<Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Water >;
};

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Voids> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Voids>
{
    using type = VoidsSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};

template<class Traits>
using TheCouplingManager = Embedded1d3dCouplingManager<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Soil> { using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::Voids>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Soil> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Soil> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Voids> { using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::Soil, TypeTag>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Voids> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Voids> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
