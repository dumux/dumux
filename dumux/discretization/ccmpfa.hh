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
 * \brief Properties for all models using cell-centered finite volume scheme with mpfa
 * \note Inherit from these properties to use a cell-centered finite volume scheme with mpfa
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_HH

#include <dune/common/reservedvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometrytraits.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumedatahandle.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cell-centered mpfa scheme.
// Create new type tags
namespace TTag {
struct CCMpfaModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the index set type used on the dual grid nodes
template<class TypeTag>
struct DualGridNodalIndexSet<TypeTag, TTag::CCMpfaModel>
{
private:
    using GV = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Traits = NodalIndexSetDefaultTraits< GV >;

public:
    using type = CCMpfaDualGridNodalIndexSet< Traits >;
};

//! Per default, we use the dynamic mpfa-o interaction volume
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::CCMpfaModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Per default, we use the dynamic mpfa-o interaction volume on boundaries
template<class TypeTag>
struct SecondaryInteractionVolume<TypeTag, TTag::CCMpfaModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::CCMpfaModel>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using PrimaryIV = GetPropType<TypeTag, Properties::PrimaryInteractionVolume>;
    using SecondaryIV = GetPropType<TypeTag, Properties::SecondaryInteractionVolume>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;
    using Traits = CCMpfaFVGridGeometryTraits<GridView, NodalIndexSet, PrimaryIV, SecondaryIV>;
public:
    using type = CCMpfaFVGridGeometry<GridView, Traits, getPropValue<TypeTag, Properties::EnableGridGeometryCache>()>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::CCMpfaModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
public:
    using type = CCMpfaGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::CCMpfaModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using FluxVariablesCacheFiller = GetPropType<TypeTag, Properties::FluxVariablesCacheFiller>;

    using PrimaryInteractionVolume = GetPropType<TypeTag, Properties::PrimaryInteractionVolume>;
    using SecondaryInteractionVolume = GetPropType<TypeTag, Properties::SecondaryInteractionVolume>;

    using PhysicsTraits = IvDataHandlePhysicsTraits<GetPropType<TypeTag, Properties::ModelTraits>>;
    using PrimaryMatVecTraits = typename PrimaryInteractionVolume::Traits::MatVecTraits;
    using SecondaryMatVecTraits = typename SecondaryInteractionVolume::Traits::MatVecTraits;

    using PrimaryIvDataHandle = InteractionVolumeDataHandle<PrimaryMatVecTraits, PhysicsTraits>;
    using SecondaryIvDataHandle = InteractionVolumeDataHandle<SecondaryMatVecTraits, PhysicsTraits>;

    using Traits = CCMpfaDefaultGridFluxVariablesCacheTraits<Problem,
                                                             FluxVariablesCache, FluxVariablesCacheFiller,
                                                             PrimaryInteractionVolume, SecondaryInteractionVolume,
                                                             PrimaryIvDataHandle, SecondaryIvDataHandle>;
public:
    using type = CCMpfaGridFluxVariablesCache<Traits, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::CCMpfaModel> { using type = CCElementBoundaryTypes; };

//! Set the BaseLocalResidual to CCLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::CCMpfaModel> { using type = CCLocalResidual<TypeTag>; };
} // namespace Properties

namespace Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethod::ccmpfa>
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
