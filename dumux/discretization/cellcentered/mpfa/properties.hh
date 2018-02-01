// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup CCMpfaDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with mpfa
 * \note Inherit from these properties to use a cell-centered finite volume scheme with mpfa
 */
#ifndef DUMUX_CC_MPFA_PROPERTIES_HH
#define DUMUX_CC_MPFA_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/mpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/connectivitymap.hh>
#include <dumux/discretization/cellcentered/mpfa/gridinteractionvolumeindexsets.hh>
#include <dumux/discretization/cellcentered/mpfa/helper.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the cell-centered mpfa scheme.
NEW_TYPE_TAG(CCMpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(CCMpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCMpfa;
};

//! Set the maximum admissible number of branches per scvf
SET_PROP(CCMpfaModel, MaxNumNeighborsPerScvf)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // Per default, we allow for 8 neighbors on network/surface grids
    static constexpr std::size_t value = dim < dimWorld ? 9 : 2;
};

//! Set the index set type used on the dual grid nodes
SET_PROP(CCMpfaModel, DualGridNodalIndexSet)
{
private:
    using GV = typename GET_PROP_TYPE(TypeTag, GridView);

    // per default, use uint8_t as iv-local index type
    using LI = std::uint8_t;

    // the specified maximum admissible number of branches per scvf
    static constexpr int maxB = GET_PROP_VALUE(TypeTag, MaxNumNeighborsPerScvf);

    // maximum admissible number of elements around a node
    // if for a given grid this number is still not high enough,
    // overwrite this property in your problem with a higher number
    static constexpr int dim = GV::dimension;
    static constexpr int maxE = dim == 3 ? 45 : 15;

public:
    using type = CCMpfaDualGridNodalIndexSet<GV, LI, dim, maxE, maxB>;
};

//! Per default, we use the dynamic mpfa-o interaction volume
SET_PROP(CCMpfaModel, PrimaryInteractionVolume)
{
public:
    using type = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
};

//! Per default, we use the mpfa-o interaction volume as secondary type
SET_PROP(CCMpfaModel, SecondaryInteractionVolume)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);

    // use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Set the default for the global finite volume geometry
SET_PROP(CCMpfaModel, FVGridGeometry)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryIV = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using SecondaryIV = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);

    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache);

    struct Traits : public DefaultMapperTraits<GridView>
    {
        using SubControlVolume = CCSubControlVolume<GridView>;
        using SubControlVolumeFace = CCMpfaSubControlVolumeFace<GridView>;
        using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);

        template< class FVGridGeometry >
        using GridIvIndexSets = CCMpfaGridInteractionVolumeIndexSets< FVGridGeometry,
                                                                      NodalIndexSet,
                                                                      PrimaryIV,
                                                                      SecondaryIV >;

        template< class FVGridGeometry, bool enableGeomCache >
        using LocalView = CCMpfaFVElementGeometry<FVGridGeometry, enableGeomCache>;

        template< class FVGridGeometry >
        using MpfaHelper = CCMpfaHelper< FVGridGeometry >;

        //! Use the correct connectivity map depending on mpfa scheme (obtain from primary iv)
        template< class FVGridGeometry >
        using ConnectivityMap = CCMpfaConnectivityMap<FVGridGeometry, FVGridGeometry::GridIVIndexSets::PrimaryInteractionVolume::MpfaMethod>;
    };
public:
    using type = CCMpfaFVGridGeometry<GridView, Traits, enableCache>;
};

//! The global flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel,
              GridFluxVariablesCache,
              CCMpfaGridFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel,
              ElementFluxVariablesCache,
              CCMpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);


//! The global previous volume variables vector class
SET_TYPE_PROP(CCMpfaModel,
              ElementVolumeVariables,
              CCMpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The global current volume variables vector class
SET_TYPE_PROP(CCMpfaModel, GridVolumeVariables, CCGridVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! Set the solution vector type for an element
SET_TYPE_PROP(CCMpfaModel, ElementSolutionVector, CCElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCMpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes<TypeTag>);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCMpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);
} // namespace Properties
} // namespace Dumux

#endif
