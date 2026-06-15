// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \brief The properties of the bump test for the bedload transport model.
 */
#ifndef DUMUX_BUMP_TEST_PROPERTIES_HH
#define DUMUX_BUMP_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/freeflow/bedload/model.hh>
#include <dumux/freeflow/shallowwater/model.hh>

#include "problem_shallowwater.hh"
#include "problem_bedload.hh"
#include "spatialparams_shallowwater.hh"
#include "spatialparams_bedload.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
    struct ShallowWaterSub{ using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
    struct BedloadSub{ using InheritsFrom = std::tuple<Bedload, CCTpfaModel>;
                       static constexpr int nGrainClasses = 1; };  // this model is limited to uniform bedload transport

} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ShallowWaterSub>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };
template<class TypeTag>
struct Grid<TypeTag, TTag::BedloadSub>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ShallowWaterSub>
{ using type = Dumux::BumpTestProblemShallowWater<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::BedloadSub>
{ using type = Dumux::BumpTestProblemBedload<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ShallowWaterSub>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = BumpTestSpatialParamsShallowWater<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BedloadSub>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = BumpTestSpatialParamsBedload<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ShallowWaterSub>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BedloadSub>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ShallowWaterSub>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::BedloadSub>
{ static constexpr bool value = false; };
} // end namespace Dumux::Properties

#endif
