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
 * \ingroup CCTpfaProperties
 * \ingroup CCTpfaModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_CCTPFA_PROPERTY_DEFAULTS_HH
#define DUMUX_CCTPFA_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/globalfvgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/globalfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/tpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>
#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

// forward declarations
template<class TypeTag> class CCElementBoundaryTypes;

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(CCTpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCTpfa;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(CCTpfaModel, GlobalFVGeometry, CCTpfaGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(CCTpfaModel, GlobalFluxVariablesCache, Dumux::CCTpfaGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(CCTpfaModel, FVElementGeometry, CCTpfaFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(CCModel, ElementVolumeVariables, Dumux::CCTpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCModel, ElementFluxVariablesCache, Dumux::CCTpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

SET_PROP(CCTpfaModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::CCTpfaSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

} // namespace Properties

} // namespace Dumux

#endif
