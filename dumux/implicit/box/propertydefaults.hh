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
 * \ingroup ImplicitProperties
 * \ingroup BoxModel
 * \file
 * \brief Default properties for box models
 */
#ifndef DUMUX_BOX_PROPERTY_DEFAULTS_HH
#define DUMUX_BOX_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/discretization/box/subcontrolvolume.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/box/globalfluxvariablescache.hh>
#include <dumux/discretization/box/elementfluxvariablescache.hh>
#include <dumux/discretization/box/globalvolumevariables.hh>
#include <dumux/discretization/box/elementvolumevariables.hh>
#include <dumux/discretization/box/globalfvgeometry.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>
#include <dumux/discretization/box/stencils.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/discretization/methods.hh>

//added for update function for fvGeometry.update in stokes/model.hh
//#include "fvelementgeometry.hh"

#include "elementboundarytypes.hh"
#include "localresidual.hh"
#include "localjacobian.hh"
#include "assembler.hh"
#include "properties.hh"

namespace Dumux {

/*!
 * \brief The box model combines the advantages of the finite-volume (FV) and finite-element (FE) methods on a dual grid
 */
// forward declarations
template<class TypeTag> class BoxLocalResidual;
template<class TypeTag> class BoxElementBoundaryTypes;
template<class TypeTag> class BoxStencilsVector;

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(BoxModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Box;
};

//! Disable evaluation of shape function gradients at the sub-control volume center by default
// The shape function gradients at the sub-control volume center are currently only
// needed for the Stokes and the linear elastic models
SET_BOOL_PROP(BoxModel, EvalGradientsAtSCVCenter, false);


//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, GlobalFVGeometry, BoxGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVElementGeometry, BoxFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The sub control volume
SET_PROP(BoxModel, SubControlVolume)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    using type = BoxSubControlVolume<ScvGeometry, IndexType>;
};

SET_PROP(BoxModel, SubControlVolumeFace)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    using type = BoxSubControlVolumeFace<ScvfGeometry, IndexType>;
};

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(BoxModel, ElementBoundaryTypes, BoxElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(BoxModel, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! The stencil container
SET_TYPE_PROP(BoxModel, StencilsVector, BoxStencilsVector<TypeTag>);

//! The global volume variables vector class
SET_TYPE_PROP(BoxModel, GlobalVolumeVariables, BoxGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The element volume variables vector class
SET_TYPE_PROP(BoxModel, ElementVolumeVariables, BoxElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(BoxModel, GlobalFluxVariablesCache, BoxGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(BoxModel, ElementFluxVariablesCache, BoxElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, BoxLocalResidual<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(BoxModel, JacobianAssembler, BoxAssembler<TypeTag>);

//! The local jacobian operator
SET_TYPE_PROP(BoxModel, LocalJacobian, BoxLocalJacobian<TypeTag>);

//! indicate that this is a box discretization
SET_BOOL_PROP(BoxModel, ImplicitIsBox, true);

} // namespace Properties
} // namespace Dumux

#endif
