// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The property specializations for the benchmark cases C1.2a/b from Schnepf et al 2023
 */
#ifndef DUMUX_1D3D_KERNEL_ROOTSOIL_PROPERTIES_HH
#define DUMUX_1D3D_KERNEL_ROOTSOIL_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/richards/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/traits.hh>

#include "problem_soil.hh"
#include "problem_root.hh"
#include "spatialparams_soil.hh"
#include "spatialparams_root.hh"

#include "couplingmanager.hh"
#include "couplingreconstruction.hh"

#ifndef BULKTYPETAG
#define BULKTYPETAG BulkCC
#endif
#ifndef LOWDIMTYPETAG
#define LOWDIMTYPETAG LowDimCC
#endif
#ifndef COUPLINGMANAGER
#define COUPLINGMANAGER CouplingManagerRootSoilKernel<Traits,CouplingReconstruction>
#endif

namespace Dumux::Properties {

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// BULK ////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
struct Bulk { using InheritsFrom = std::tuple<Richards>; };
struct BulkCC { using InheritsFrom = std::tuple<Bulk, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Bulk> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Bulk> { using type = SoilProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Bulk>
{ using type = SoilSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::Bulk> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::Bulk> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::Bulk> { static constexpr bool value = true; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Bulk> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::Bulk> { static constexpr bool value = false; };

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// EMBEDDED ////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
struct LowDim { using InheritsFrom = std::tuple<OneP>; };
struct LowDimCC { using InheritsFrom = std::tuple<LowDim, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LowDim> { using type = Dune::FoamGrid<1, 3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LowDim> { using type = RootProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LowDim>
{ using type = FluidSystems::OnePLiquid<double, Components::SimpleH2O<double> >; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LowDim>
{ using type = RootSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

template<class TypeTag> struct EnableGridGeometryCache<TypeTag, TTag::LowDim> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridVolumeVariablesCache<TypeTag, TTag::LowDim> { static constexpr bool value = true; };
template<class TypeTag> struct EnableGridFluxVariablesCache<TypeTag, TTag::LowDim> { static constexpr bool value = true; };
template<class TypeTag> struct SolutionDependentAdvection<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentMolecularDiffusion<TypeTag, TTag::LowDim> { static constexpr bool value = false; };
template<class TypeTag> struct SolutionDependentHeatConduction<TypeTag, TTag::LowDim> { static constexpr bool value = false; };

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// COUPLING ////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

template<class Traits, class CouplingReconstruction>
using TheCouplingManager = COUPLINGMANAGER;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BULKTYPETAG>
{
    using Traits = MultiDomainTraits<TTag::BULKTYPETAG, TTag::LOWDIMTYPETAG>;
    using CouplingReconstruction = Dumux::RootSoil::CouplingReconstruction<typename GetPropType<TypeTag, Properties::SpatialParams>::PcKrSwCurve>;
    using type = TheCouplingManager<Traits, CouplingReconstruction>;
};
template<class TypeTag>
struct PointSource<TypeTag, TTag::BULKTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::BULKTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LOWDIMTYPETAG>
{
    using Traits = MultiDomainTraits<TTag::BULKTYPETAG, TTag::LOWDIMTYPETAG>;
    using CouplingReconstruction = Dumux::RootSoil::CouplingReconstruction<typename GetPropType<Properties::TTag::BULKTYPETAG, Properties::SpatialParams>::PcKrSwCurve>;
    using type = TheCouplingManager<Traits, CouplingReconstruction>;
};
template<class TypeTag>
struct PointSource<TypeTag, TTag::LOWDIMTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::LOWDIMTYPETAG>
{ using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
