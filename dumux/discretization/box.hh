// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_DISCRETIZTAION_BOX_HH
#define DUMUX_DISCRETIZTAION_BOX_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/discretization/localdoftraits.hh>

#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvolumevariables.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux::Properties {

//! Type tag for the box scheme.
// Create new type tags
namespace TTag {
struct BoxModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::BoxModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::BoxModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Traits = CVFEDefaultGridVolumeVariablesTraits<Problem, VolumeVariables>;
public:
    using type = CVFEGridVolumeVariables<Traits, enableCache>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::BoxModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
public:
    using type = CVFEGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache>;
};

//! The flux variables cache type
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::BoxModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CVFEFluxVariablesCache<Scalar, GridGeometry>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::BoxModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    using type = CVFEElementBoundaryTypes<BoundaryTypes>;
};

//! Set the BaseLocalResidual to CVFELocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::BoxModel>
{ using type = CVFELocalResidual<TypeTag>; };

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::Box>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GG::SubControlVolume;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scv)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolume>()))>;
};

template<class GridView>
struct LocalDofTraits<GridView, DiscretizationMethods::Box>
{
    static constexpr int dim = GridView::dimension;
    // Dofs are located at the corners
    static constexpr int numCubeElementDofs = (1<<dim);
};

} // end namespace Dumux:Detail

#endif
