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
 * \ingroup CCMpfaProperties
 * \ingroup CCMpfaModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_CCMPFA_PROPERTY_DEFAULTS_HH
#define DUMUX_CCMPFA_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/globalfvgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/globalvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/globalfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/mpfa/helper.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolume.hh>
#include <dumux/discretization/cellcentered/mpfa/globalinteractionvolumeseeds.hh>
#include <dumux/discretization/cellcentered/mpfa/interiorboundarydata.hh>
#include <dumux/implicit/cellcentered/mpfa/localresidual.hh>
#include <dumux/implicit/cellcentered/properties.hh>

namespace Dumux {

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(CCMpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCMpfa;
};

//! Set the BaseLocalResidual to CCMpfaLocalResidual
SET_TYPE_PROP(CCMpfaModel, BaseLocalResidual, CCMpfaLocalResidual<TypeTag>);

//! By default we set the o-method as the Mpfa method of choice
SET_PROP(CCMpfaModel, MpfaMethod)
{
    static const MpfaMethods value = MpfaMethods::oMethod;
};

//! The mpfa helper class
SET_TYPE_PROP(CCMpfaModel, MpfaHelper, CCMpfaHelper<TypeTag>);

//! The interaction volume class
SET_TYPE_PROP(CCMpfaModel, InteractionVolume, CCMpfaInteractionVolume<TypeTag>);

//! The boundary interaction volume class (for methods other than the omethod)
SET_TYPE_PROP(CCMpfaModel, BoundaryInteractionVolume, typename GET_PROP_TYPE(TypeTag, InteractionVolume)::BoundaryInteractionVolume);

//! The global interaction volume seeds class
SET_TYPE_PROP(CCMpfaModel, GlobalInteractionVolumeSeeds, CCMpfaGlobalInteractionVolumeSeeds<TypeTag>);

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(CCMpfaModel, GlobalFVGeometry, CCMpfaGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The global current volume variables vector class
SET_TYPE_PROP(CCMpfaModel, GlobalVolumeVariables, Dumux::CCMpfaGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel, GlobalFluxVariablesCache, CCMpfaGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(CCMpfaModel, FVElementGeometry, CCMpfaFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(CCMpfaModel, ElementVolumeVariables, Dumux::CCMpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel, ElementFluxVariablesCache, Dumux::CCMpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

SET_PROP(CCMpfaModel, SubControlVolumeFace)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    typedef Dumux::CCMpfaSubControlVolumeFace<GET_PROP_VALUE(TypeTag, MpfaMethod), ScvfGeometry, IndexType> type;
};

// By default, we use tpfa on the boundaries
SET_BOOL_PROP(CCMpfaModel, UseTpfaBoundary, true);

// By default, we disable interior boundaries
SET_BOOL_PROP(CCMpfaModel, EnableInteriorBoundaries, false);

// By default, interior boundaries are static
SET_BOOL_PROP(CCMpfaModel, MpfaFacetCoupling, false);

// The default interior Dirichlet boundary data
SET_TYPE_PROP(CCMpfaModel, InteriorBoundaryData, InteriorBoundaryData<TypeTag>);

// By default, we use simple coupling conditions (Xi = 1)
SET_SCALAR_PROP(CCMpfaModel, MpfaXi, 1.0);

// By default, we set the quadrature point to the mid point of the element facets
SET_SCALAR_PROP(CCMpfaModel, MpfaQ, 0.0);

// By default, the auxiliary volume size is half of the interaction region
SET_SCALAR_PROP(CCMpfaModel, MpfaC, 0.5);

// By default, we set the quadrature point 1/2 away from the mid point of the auxiliary sub faces
SET_SCALAR_PROP(CCMpfaModel, MpfaP, 0.75);

} // namespace Properties

} // namespace Dumux

#endif
