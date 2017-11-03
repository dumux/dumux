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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tag and some properties for models using
 *        a cell-centered scheme with two-point flux approximation.
 */

#ifndef DUMUX_CC_TPFA_PROPERTIES_HH
#define DUMUX_CC_TPFA_PROPERTIES_HH

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/globalvolumevariables.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>

#include <dumux/implicit/cellcentered/elementboundarytypes.hh>
#include <dumux/implicit/cellcentered/localresidual.hh>

#include <dumux/discretization/cellcentered/connectivitymap.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/globalfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/tpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the cell-centered tpfa scheme.
NEW_TYPE_TAG(CCTpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(CCTpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCTpfa;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(CCTpfaModel, FVGridGeometry, CCTpfaFVGridGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(CCTpfaModel, GlobalFluxVariablesCache, CCTpfaGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(CCTpfaModel, FVElementGeometry, CCTpfaFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(CCTpfaModel, ElementVolumeVariables, CCTpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCTpfaModel, ElementFluxVariablesCache, CCTpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The global current volume variables vector class
SET_TYPE_PROP(CCTpfaModel, GlobalVolumeVariables, CCGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The sub control volume
SET_PROP(CCTpfaModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    using type = CCSubControlVolume<ScvGeometry, IndexType>;
};

//! The sub control volume face
SET_PROP(CCTpfaModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::CCTpfaSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! Set the solution vector type for an element
SET_TYPE_PROP(CCTpfaModel, ElementSolutionVector, CCElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCTpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes<TypeTag>);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCTpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);

//! Set the AssemblyMap to the SimpleAssemblyMap per default
SET_TYPE_PROP(CCTpfaModel, AssemblyMap, CCSimpleConnectivityMap<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
