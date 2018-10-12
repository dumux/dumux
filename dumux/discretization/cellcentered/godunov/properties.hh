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
 * \ingroup GodunovDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with Godunov
 * \note Inherit from these properties to use a cell-centered finite volume scheme with Godunov
 */

#ifndef DUMUX_CC_Godunov_PROPERTIES_HH
#define DUMUX_CC_Godunov_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/discretization/cellcentered/godunov/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/godunov/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/godunov/subcontrolvolumeface.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cell-centered godunov scheme.
NEW_TYPE_TAG(GodunovModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the default for the global finite volume geometry
SET_PROP(GodunovModel, FVGridGeometry)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = GodunovFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
SET_PROP(GodunovModel, GridVolumeVariables)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    using type = GodunovGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};


//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(GodunovModel, ElementBoundaryTypes, CCElementBoundaryTypes);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(GodunovModel, BaseLocalResidual, CCLocalResidual<TypeTag>);
} // namespace Properties
} // namespace Dumux

#endif
