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
#ifndef DUMUX_STAGGERED_PROPERTIES_HH
#define DUMUX_STAGGERED_PROPERTIES_HH

#include <dumux/implicit/properties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup StaggeredModel
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the StaggeredModel.
 */
namespace Dumux
{

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on staggered scheme
NEW_TYPE_TAG(StaggeredModel, INHERITS_FROM(ImplicitBase));

NEW_PROP_TAG(StaggeredGeometryHelper); //!< Helper class to ease the creation of stencils, etc.

NEW_PROP_TAG(CellCenterPrimaryVariables); //!< A vector of primary variables for cell center dofs
NEW_PROP_TAG(FacePrimaryVariables); //!< A vector of primary variables for face dofs
NEW_PROP_TAG(CellCenterSolutionVector); //!< Vector containing all cell centered primary variables
NEW_PROP_TAG(FaceSolutionVector); //!< Vector containing all face primary variables

NEW_PROP_TAG(EnableInteriorBoundaries); //!< For compatibility

NEW_PROP_TAG(BaseEpsilon); //!< Set one or different base epsilons for the calculations of the localJacobian's derivatives

NEW_PROP_TAG(NumEqCellCenter); //!< Number of equations per cell center dof
NEW_PROP_TAG(NumEqFace); //!< Number of equations per face dof
NEW_PROP_TAG(DofTypeIndices); //!< Indices to choose between cell center and face dofs

NEW_PROP_TAG(FaceVariables); //!< Variables associated to facets (equivalent to volVars)
NEW_PROP_TAG(GlobalFaceVars); //!< Global vector of facet variables

NEW_PROP_TAG(VtkWriteFaceData); //!< Decide whether to write separate vtp files for face variables

}
}

// \}

#include <dumux/implicit/staggered/propertydefaults.hh>

#endif
