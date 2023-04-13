// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>
#if HAVE_QUADMATH
#include <dune/common/quadmath.hh>
#endif
#include <dumux/common/boundarytypes.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1bubble.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2>
#endif

////////////////////////////////////////////////////
// A fake grid that cannot communicate.           //
// Can be used to make sure that compilation and  //
// sequential execution also work for such grids. //
////////////////////////////////////////////////////
namespace Dumux {

template<int dim>
class NoCommunicateGrid;

template<int dim>
class NoCommunicateGridLeafGridView
: public Dune::YaspGrid<dim>::LeafGridView
{
    using ParentType = typename Dune::YaspGrid<dim>::LeafGridView;
public:
    using ParentType::ParentType;

    struct Traits : public ParentType::Traits
    { using Grid = NoCommunicateGrid<dim>; };

    const typename Traits::Grid& grid() const
    { return static_cast<const typename Traits::Grid&>(ParentType::grid()); }
};

template<int dim>
class NoCommunicateGrid : public Dune::YaspGrid<dim>
{
    using ParentType = Dune::YaspGrid<dim>;
public:
    using ParentType::ParentType;
    struct Traits : public ParentType::Traits
    { using LeafGridView = NoCommunicateGridLeafGridView<dim>; };

    using LeafGridView = NoCommunicateGridLeafGridView<dim>;

    typename Traits::LeafGridView leafGridView() const
    { return NoCommunicateGridLeafGridView<dim>(*this); }
private:
    using ParentType::communicate;
    using ParentType::communicateCodim;
};

template<int dim>
class GridManager<NoCommunicateGrid<dim>> : public GridManager<Dune::YaspGrid<dim>>
{
    using ParentType = GridManager<Dune::YaspGrid<dim>>;
public:
    using ParentType::ParentType;
    using Grid = NoCommunicateGrid<dim>;
    Grid& grid() { return static_cast<NoCommunicateGrid<dim>&>(ParentType::grid()); }
};

} // end namespace Dumux

namespace Dune::Capabilities {

template<int dim, int codim>
struct canCommunicate<Dumux::NoCommunicateGrid<dim>, codim>
{ static constexpr bool v = false; };

} // end namespace Dune::Capabilities
/////////////////////////////////////////

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressibleTpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCTpfaModel>; };
struct OnePIncompressibleMpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCMpfaModel>; };
struct OnePIncompressibleBox { using InheritsFrom = std::tuple<OnePIncompressible, BoxModel>; };
struct OnePIncompressibleFCDiamond { using InheritsFrom = std::tuple<OnePIncompressible, FaceCenteredDiamondModel>; };
struct OnePIncompressiblePQ1Bubble { using InheritsFrom = std::tuple<OnePIncompressible, PQ1BubbleModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressible> { using type = GRIDTYPE; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressible>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePIncompressible> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };

// define a TypeTag for a quad precision test
#if HAVE_QUADMATH
namespace TTag {
struct OnePIncompressibleTpfaQuad { using InheritsFrom = std::tuple<OnePIncompressibleTpfa>; };
} // end namespace TTag
template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePIncompressibleTpfaQuad> { using type = Dune::Float128; };
#endif
}
#endif
