// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The properties for the single phase + tracer test
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_PROPERTIES_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_PROPERTIES_HH

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/tracer/model.hh>

#include "spatialparams_1p.hh"
#include "spatialparams_tracer.hh"
#include "problem_1p_bulk.hh"
#include "problem_1p_lowdim.hh"
#include "problem_tracer_bulk.hh"
#include "problem_tracer_lowdim.hh"
#include "tracerfluidsystem.hh"
#include "tracermodeltraits.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::SimpleH2O<Scalar> >;
};

// Create new type tags
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimMpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimBox { using InheritsFrom = std::tuple<OnePLowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = Dune::FoamGrid<1, 2>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::SimpleH2O<Scalar> >;
};

// Create new type tags
namespace TTag {
struct TracerTestBulk { using InheritsFrom = std::tuple<Tracer>; };

// define the type tags
struct TracerBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, TracerTestBulk>; };
struct TracerBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, TracerTestBulk>; };
struct TracerBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, TracerTestBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTestBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

//! Overwrite the advection type property
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::TracerBulkTpfa> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::TracerBulkBox> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTestBulk> { using type = TracerBulkProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTestBulk>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTestBulk> { static constexpr bool value = false; };

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::TracerTestBulk> { static constexpr bool value = false; };

//! set the model traits (with disabled diffusion)
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TracerTestBulk>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = TracerTestModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(),
                                                                   getPropValue<TypeTag, Properties::EnableCompositionalDispersion>()>;
};

// use the test-specific fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTestBulk> { using type = TracerFluidSystem<TypeTag>; };

// Create new type tags
namespace TTag {
struct TracerTestLowDim { using InheritsFrom = std::tuple<Tracer>; };

// define the type tags for both bulk and lowdim type tag here
struct TracerLowDimTpfa { using InheritsFrom = std::tuple<TracerTestLowDim, CCTpfaModel>; };
struct TracerLowDimMpfa { using InheritsFrom = std::tuple<TracerTestLowDim, CCTpfaModel>; };
struct TracerLowDimBox { using InheritsFrom = std::tuple<TracerTestLowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTestLowDim> { using type = Dune::FoamGrid<1, 2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTestLowDim> { using type = TracerLowDimProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTestLowDim>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTestLowDim> { static constexpr bool value = false; };

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::TracerTestLowDim> { static constexpr bool value = false; };

//! set the model traits (with disabled diffusion)
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TracerTestLowDim>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = TracerTestModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(),
                                                                   getPropValue<TypeTag, Properties::EnableCompositionalDispersion>()>;
};

// use the test-specific fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTestLowDim> { using type = TracerFluidSystem<TypeTag>; };

template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = Dumux::GetPropType<BulkTypeTag, Dumux::Properties::GridGeometry>;
    using LowDimFVGridGeometry = Dumux::GetPropType<LowDimTypeTag, Dumux::Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property for the box test
using BoxTraits = TestTraits<Properties::TTag::OnePBulkBox, Properties::TTag::OnePLowDimBox>;
using BoxTracerTraits = TestTraits<Properties::TTag::TracerBulkBox, Properties::TTag::TracerLowDimBox>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkBox> { using type = typename BoxTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimBox> { using type = typename BoxTracerTraits::CouplingManager; };

// set cm property for the tpfa test
using TpfaTraits = TestTraits<Properties::TTag::OnePBulkTpfa, Properties::TTag::OnePLowDimTpfa>;
using TpfaTracerTraits = TestTraits<Properties::TTag::TracerBulkTpfa, Properties::TTag::TracerLowDimTpfa>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkTpfa> { using type = typename TpfaTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimTpfa> { using type = typename TpfaTracerTraits::CouplingManager; };

// set cm property for the mpfa test
using MpfaTraits = TestTraits<Properties::TTag::OnePBulkMpfa, Properties::TTag::OnePLowDimMpfa>;
using MpfaTracerTraits = TestTraits<Properties::TTag::TracerBulkMpfa, Properties::TTag::TracerLowDimMpfa>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkMpfa> { using type = typename MpfaTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimMpfa> { using type = typename MpfaTracerTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
