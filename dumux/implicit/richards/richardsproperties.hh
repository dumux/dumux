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
 * \ingroup RichardsModel
 * \file
 *
 * \brief Contains the property declarations for the Richards box
 *        model.
 */
#ifndef DUMUX_RICHARDS_PROPERTIES_HH
#define DUMUX_RICHARDS_PROPERTIES_HH

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/implicit/nonisothermal/niproperties.hh>

namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal richards model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
NEW_TYPE_TAG(Richards);
NEW_TYPE_TAG(BoxRichards, INHERITS_FROM(BoxModel, Richards));
NEW_TYPE_TAG(CCRichards, INHERITS_FROM(CCModel, Richards));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(RichardsNI, INHERITS_FROM(Richards, NonIsothermal));
NEW_TYPE_TAG(BoxRichardsNI, INHERITS_FROM(BoxModel, RichardsNI));
NEW_TYPE_TAG(CCRichardsNI, INHERITS_FROM(CCModel, RichardsNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (by default extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The type of the parameter object for the material law (by default extracted from the spatial parameters)
NEW_PROP_TAG(FluidSystem); //!< The fluid system to be used for the Richards model
NEW_PROP_TAG(WettingPhase); //!< Fluid which represents the wetting phase
NEW_PROP_TAG(NonwettingPhase); //!< Fluid which represents the non-wetting phase
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the weight of the upwind direction in the mass conservation equations
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< The value of the weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(UseHead); //!< Defines whether pressure [Pa] (false) or pressure head [cm] (true) is used

// \}
}

} // end namespace

#endif
