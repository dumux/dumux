// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the staggered scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 * TODO: detailed documentation and figures
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_HH

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/fclocalresidual.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/staggered/gridvolumevariables.hh>
#include <dumux/discretization/facecentered/staggered/gridfluxvariablescache.hh>
#include <dumux/discretization/facecentered/staggered/elementboundarytypes.hh>

namespace Dumux::Properties {

//! Type tag for the staggered scheme.
// Create new type tags
namespace TTag {
struct FaceCenteredStaggeredModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
public:
    using type = FaceCenteredStaggeredFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Traits = FaceCenteredStaggeredDefaultGridVolumeVariablesTraits<Problem, VolumeVariables>;
public:
    using type = FaceCenteredStaggeredGridVolumeVariables<Traits, enableCache>;
};


//! Set the global flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
    using FluxVariablesCacheFiller = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCacheFiller, FluxVariablesCaching::EmptyCacheFiller
    >;
public:
    using type = FaceCenteredStaggeredGridFluxVariablesCache<Problem, FluxVariablesCache, FluxVariablesCacheFiller, enableCache>;
};

//! Set the grid variables (volume, flux and face variables)
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFVC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
public:
    using type = FVGridVariables<GG, GVV, GFVC>;
};

//! Set the BaseLocalResidual to FaceCenteredLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::FaceCenteredStaggeredModel> { using type = FaceCenteredLocalResidual<TypeTag>; };

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    using type = FaceCenteredStaggeredElementBoundaryTypes<BoundaryTypes>;
};

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::FCStaggered>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scv)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolumeFace>()))>;
};

} // end namespace Detail

#endif
