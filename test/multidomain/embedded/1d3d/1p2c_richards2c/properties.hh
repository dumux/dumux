// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The properties for the one-phase soil problem.
 */
#ifndef DUMUX_TISSUE_PROPERTIES_HH
#define DUMUX_TISSUE_PROPERTIES_HH



#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/richardsnc/model.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_soil.hh"
#include "problem_root.hh"

#include "spatialparams_soil.hh"
#include "spatialparams_root.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Soil { using InheritsFrom = std::tuple<RichardsNC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::Soil> { using type = Dune::UGGrid<3>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Soil> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>; };
#endif

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Soil> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Soil> { using type = SoilProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Soil>
{
    using type = SoilSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                   GetPropType<TypeTag, Properties::Scalar>>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Soil>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Soil> { static constexpr bool value = true; };

// Create new type tags
namespace TTag {
struct Root { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Root> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Root> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Root> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Root> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Root> { using type = RootProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Root>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Root>
{
    using type = RootSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                   GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Root> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Soil>
{
    using Traits = MultiDomainTraits<TypeTag, Properties::TTag::Root>;
    using type = Embedded1d3dCouplingManager<Traits, Embedded1d3dCouplingMode::Average>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Root>
{
    using Traits = MultiDomainTraits<Properties::TTag::Soil, TypeTag>;
    using type = Embedded1d3dCouplingManager<Traits, Embedded1d3dCouplingMode::Average>;
};

template<class TypeTag>
struct PointSource<TypeTag, TTag::Soil> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Root> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Soil> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Root> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
