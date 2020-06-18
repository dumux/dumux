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
 * \brief Defines a type tag and some properties for models using the staggered scheme.
          This scheme features degrees of freedom at the elements' centers and intersections (faces).
 * TODO: detailed documentation and figures
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERD_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERD_HH

#include <dumux/common/properties.hh>

#include <dumux/assembly/fclocalresidual.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/staggered/gridvolumevariables.hh>
#include <dumux/discretization/facecentered/staggered/gridfluxvariablescache.hh>
#include <dumux/discretization/facecentered/staggered/elementboundarytypes.hh>



namespace Dumux {

// forward declarations
class CCElementBoundaryTypes;

namespace Properties
{
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
// TODO
//! Set the grid variables (volume, flux and face variables)
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::FaceCenteredStaggeredModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
public:
    using type = FaceCenteredStaggeredGridFluxVariablesCache<Problem, FluxVariablesCache>;
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

//! Use the cell center element boundary types per default
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::FaceCenteredStaggeredModel> { using type = FaceCenteredStaggeredElementBoundaryTypes<GetPropType<TypeTag, Properties::BoundaryTypes>>; };

//! Boundary types at a single degree of freedom
template<class TypeTag>
struct BoundaryTypes<TypeTag, TTag::FaceCenteredStaggeredModel> { using type = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

//! Set the BaseLocalResidual to BoxLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::FaceCenteredStaggeredModel> { using type = FaceCenteredLocalResidual<TypeTag>; };

// TODO: bundle SolutionVector, JacobianMatrix
//       in LinearAlgebra traits

} // namespace Properties
} // namespace Dumux

#endif
