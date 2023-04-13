// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The properties for a fracture problem.
 */
#ifndef DUMUX_FRACTURE_PROPERTIES_HH
#define DUMUX_FRACTURE_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/embedded/couplingmanager2d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "spatialparams.hh"
#include "problem_matrix.hh"
#include "problem_fracture.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Matrix { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Matrix> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 3> >; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Matrix> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Matrix> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Matrix> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Matrix> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Matrix> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Matrix> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Matrix> { using type = MatrixProblem<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Matrix> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Matrix>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Matrix>
{
    using type = MatrixFractureSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                             GetPropType<TypeTag, Properties::Scalar>>;
};

// Create new type tags
namespace TTag {
struct Fracture { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Fracture> { using type = Dune::FoamGrid<2, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Fracture> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Fracture> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Fracture> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Fracture> { using type = FractureProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Fracture>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Fracture> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Fracture>
{
    using type = MatrixFractureSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                             GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Matrix>
{
    using Traits = MultiDomainTraits<TypeTag, Properties::TTag::Fracture>;
    using type = EmbeddedCouplingManager2d3d<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::Fracture>
{
    using Traits = MultiDomainTraits<Properties::TTag::Matrix, TypeTag>;
    using type = EmbeddedCouplingManager2d3d<Traits>;
};

template<class TypeTag>
struct PointSource<TypeTag, TTag::Matrix> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::Fracture> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Matrix> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::Fracture> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
