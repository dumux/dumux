// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
 * \brief Default properties for box models
 */
#ifndef DUMUX_CC_PROPERTY_DEFAULTS_HH
#define DUMUX_CC_PROPERTY_DEFAULTS_HH

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dumux/linear/ovlpistlsolverbackend.hh>
#include <dumux/linear/amgbackend.hh>
#endif

#include <dumux/implicit/common/implicitpropertydefaults.hh>
#include "ccassembler.hh"
#include "ccfvelementgeometry.hh"
#include "ccelementboundarytypes.hh"
#include "cclocalresidual.hh"
#include "ccelementvolumevariables.hh"
#include "ccproperties.hh"

namespace Dumux {

// forward declarations
template<class TypeTag> class CCModel;
template<class TypeTag> class CCLocalResidual;
template<class TypeTag> class CCElementBoundaryTypes;
template<class TypeTag> class CCElementVolumeVariables;
template<class TypeTag> class CCFVElementGeometry;

namespace Properties {
//! Set the default for the FVElementGeometry
SET_TYPE_PROP(CCModel, FVElementGeometry, Dumux::CCFVElementGeometry<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(CCModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCModel, BaseLocalResidual, Dumux::CCLocalResidual<TypeTag>);

//! An array of secondary variable containers
SET_TYPE_PROP(CCModel, ElementVolumeVariables, Dumux::CCElementVolumeVariables<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(CCModel, JacobianAssembler, Dumux::CCAssembler<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(CCModel, ImplicitIsBox, false);

#if HAVE_DUNE_PDELAB
//! use the element-wise constant local FEM space by default
SET_PROP(CCModel, ImplicitLocalFemSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::P0LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

SET_PROP(CCModel, ImplicitPDELabBackend)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef Dumux::ISTLBackend_BCGS_AMG_SSOR<GridOperator> type;
};
#endif

} // namespace Properties
} // namespace Dumux

#endif
