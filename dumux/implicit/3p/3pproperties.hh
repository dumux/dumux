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
 * \ingroup ThreePModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the three-phase fully implicit model.
 */
#ifndef DUMUX_3P_PROPERTIES_HH
#define DUMUX_3P_PROPERTIES_HH

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/implicit/nonisothermal/niproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tags for the implicit isothermal one-phase two-component problems
NEW_TYPE_TAG(ThreeP);
NEW_TYPE_TAG(BoxThreeP, INHERITS_FROM(BoxModel, ThreeP));
NEW_TYPE_TAG(CCThreeP, INHERITS_FROM(CCModel, ThreeP));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(ThreePNI, INHERITS_FROM(ThreeP, NonIsothermal));
NEW_TYPE_TAG(BoxThreePNI, INHERITS_FROM(BoxModel, ThreePNI));
NEW_TYPE_TAG(CCThreePNI, INHERITS_FROM(CCModel, ThreePNI));


//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(FluidState); //!< the phases' state

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The parameters of the material law (extracted from the spatial parameters)

NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(BaseFluxVariables); //! The base flux variables
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}
}

#endif
