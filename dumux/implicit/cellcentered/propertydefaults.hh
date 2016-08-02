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
 * \ingroup CCProperties
 * \ingroup CCModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_CC_PROPERTY_DEFAULTS_HH
#define DUMUX_CC_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/discretization/cellcentered/globalvolumevariables.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/stencils.hh>

#include "elementboundarytypes.hh"
#include "localresidual.hh"
#include "properties.hh"
#include "localjacobian.hh"
#include "assembler.hh"

namespace Dumux
{

// forward declaration
template<class TypeTag> class CCElementBoundaryTypes;
template<class TypeTag> class CCLocalResidual;
template<class TypeTag> class CCStencilsVector;

namespace Properties
{
//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCModel, ElementBoundaryTypes, CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(CCModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The local jacobian operator
SET_TYPE_PROP(CCModel, LocalJacobian, Dumux::CCLocalJacobian<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(CCModel, JacobianAssembler, Dumux::CCAssembler<TypeTag>);

//! The stencil container
SET_TYPE_PROP(CCModel, StencilsVector, CCStencilsVector<TypeTag>);

//! The sub control volume
SET_PROP(CCModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::CCSubControlVolume<ScvGeometry, IndexType> type;
};

//! The global current volume variables vector class
SET_TYPE_PROP(CCModel, GlobalVolumeVariables, Dumux::CCGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCModel, BaseLocalResidual, CCLocalResidual<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(CCModel, ImplicitIsBox, false);

} // end namespace Properties
} // end namespace Dumux

#endif
