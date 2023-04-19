// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Properties for all models using cell-centered finite volume scheme with TPFA
 * \note Inherit from these properties to use a cell-centered finite volume scheme with TPFA
 */

#ifndef DUMUX_DISCRETIZATION_CC_TPFA_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/tpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cell-centered tpfa scheme.
// Create new type tags
namespace TTag {
struct CCTpfaModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::CCTpfaModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::CCTpfaModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
public:
    using type = CCTpfaGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::CCTpfaModel>
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
    using type = CCTpfaGridFluxVariablesCache<Problem, FluxVariablesCache, FluxVariablesCacheFiller, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::CCTpfaModel> { using type = CCElementBoundaryTypes; };

//! Set the BaseLocalResidual to CCLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::CCTpfaModel> { using type = CCLocalResidual<TypeTag>; };
} // namespace Properties

namespace Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::CCTpfa>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scvf)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolumeFace>()))>;
};

} // end namespace Detail

} // namespace Dumux

#endif
