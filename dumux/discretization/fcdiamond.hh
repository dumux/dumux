// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the diamond scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_HH

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/fcdiamondlocalresidual.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/diamond/gridvolumevariables.hh>
#include <dumux/discretization/facecentered/diamond/gridfluxvariablescache.hh>
#include <dumux/discretization/facecentered/diamond/fluxvariablescache.hh>
#include <dumux/discretization/facecentered/diamond/elementboundarytypes.hh>

namespace Dumux::Properties {

//! Type tag for the staggered scheme.
// Create new type tags
namespace TTag {
struct FaceCenteredDiamondModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
public:
    using type = FaceCenteredDiamondFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Traits = FaceCenteredDiamondDefaultGridVolumeVariablesTraits<Problem, VolumeVariables>;
public:
    using type = FaceCenteredDiamondGridVolumeVariables<Traits, enableCache>;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FaceCenteredDiamondFluxVariablesCache<Scalar, GridGeometry>;
};

//! Set the global flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
    using FluxVariablesCacheFiller = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCacheFiller, FluxVariablesCaching::EmptyCacheFiller
    >;
public:
    using type = FaceCenteredDiamondGridFluxVariablesCache<Problem, GridGeometry, FluxVariablesCache, FluxVariablesCacheFiller, enableCache>;
};

//! Set the grid variables (volume, flux and face variables)
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GVV = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using GFVC = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
public:
    using type = FVGridVariables<GG, GVV, GFVC>;
};

//! Set the BaseLocalResidual to DiamondLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::FaceCenteredDiamondModel>
{ using type = FaceCenteredDiamondLocalResidual<TypeTag>; };

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::FaceCenteredDiamondModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
    using FVGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
public:
    using type = FaceCenteredDiamondElementBoundaryTypes<BoundaryTypes, FVGeometry>;
};

} // namespace Dumux::Properties

namespace Dumux::Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::FCDiamond>
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

} // end namespace Dumux::Detail

#endif
