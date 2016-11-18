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
 * \ingroup StaggeredModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_STAGGERED_PROPERTY_DEFAULTS_HH
#define DUMUX_STAGGERED_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
// #include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/discretization/staggered/globalfvgeometry.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>
#include <dumux/implicit/staggered/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/staggered/stencils.hh>


#include <dumux/freeflow/staggered/fluxvariables.hh>
#include <dumux/freeflow/staggered/fluxvariablescache.hh>




#include "assembler.hh"
#include "localresidual.hh"
#include "localjacobian.hh"

namespace Dumux {

// forward declarations
template<class TypeTag> class CCElementBoundaryTypes;

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(StaggeredModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Staggered;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(StaggeredModel, GlobalFVGeometry, StaggeredGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(StaggeredModel, FVElementGeometry, StaggeredFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The sub control volume
SET_PROP(StaggeredModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::StaggeredSubControlVolume<ScvGeometry, IndexType> type;
};

SET_PROP(StaggeredModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::StaggeredSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(StaggeredModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(StaggeredModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The global current volume variables vector class
SET_TYPE_PROP(StaggeredModel, GlobalVolumeVariables, Dumux::CCGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(StaggeredModel, GlobalFluxVariablesCache, Dumux::CCGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The local jacobian operator
SET_TYPE_PROP(StaggeredModel, LocalJacobian, Dumux::StaggeredLocalJacobian<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(StaggeredModel, JacobianAssembler, Dumux::StaggeredAssembler<TypeTag>);

//! The stencil container
SET_TYPE_PROP(StaggeredModel, StencilsVector, Dumux::StaggeredStencilsVector<TypeTag>);

//! The local flux variables cache vector class
SET_TYPE_PROP(StaggeredModel, ElementFluxVariablesCache, Dumux::CCElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(StaggeredModel, ElementVolumeVariables, Dumux::CCElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! Set the BaseLocalResidual to StaggeredLocalResidual
SET_TYPE_PROP(StaggeredModel, BaseLocalResidual, Dumux::StaggeredLocalResidual<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(StaggeredModel, ImplicitIsBox, false);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(StaggeredModel, FluxVariables, FreeFlowFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(StaggeredModel, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

} // namespace Properties

} // namespace Dumux

#endif
