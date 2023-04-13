// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The properties for the Poiseuille flow problem.
 */
#ifndef DUMUX_POISEUILLE_FLOW_TEST_PROPERTIES_HH
#define DUMUX_POISEUILLE_FLOW_TEST_PROPERTIES_HH

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>
#endif

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

/*!
 * \ingroup ShallowWaterTests
 * \brief The properties for the Poiseuille flow problem test.
 *
 */
namespace TTag {
struct PoiseuilleFlow { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::PoiseuilleFlow> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoiseuilleFlow> { using type = Dumux::PoiseuilleFlowProblem<TypeTag> ; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoiseuilleFlow>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

public:
    using type = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoiseuilleFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PoiseuilleFlow> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
