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
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_BOX_PROPERTIES_HH
#define DUMUX_BOX_PROPERTIES_HH

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/implicit/box/elementboundarytypes.hh>
#include <dumux/implicit/box/localresidual.hh>

#include <dumux/discretization/box/subcontrolvolume.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/box/globalfluxvariablescache.hh>
#include <dumux/discretization/box/elementfluxvariablescache.hh>
#include <dumux/discretization/box/globalvolumevariables.hh>
#include <dumux/discretization/box/elementvolumevariables.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the box scheme.
NEW_TYPE_TAG(BoxModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(BoxModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Box;
};

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVGridGeometry, BoxFVGridGeometry<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVElementGeometry, BoxFVElementGeometry<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

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

//! Set the solution vector type for an element
SET_TYPE_PROP(BoxModel, ElementSolutionVector, BoxElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(BoxModel, ElementBoundaryTypes, BoxElementBoundaryTypes<TypeTag>);

//! The global volume variables vector class
SET_TYPE_PROP(BoxModel, GlobalVolumeVariables, BoxGlobalVolumeVariables<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The element volume variables vector class
SET_TYPE_PROP(BoxModel, ElementVolumeVariables, BoxElementVolumeVariables<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(BoxModel, GlobalFluxVariablesCache, BoxGlobalFluxVariablesCache<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(BoxModel, ElementFluxVariablesCache, BoxElementFluxVariablesCache<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, BoxLocalResidual<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
