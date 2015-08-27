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
 * \ingroup TwoPNCModel
 *
 * \file
 *
 * \brief Defines the properties required for the two-phase n-component
 *        fully implicit model.
 */
#ifndef DUMUX_2PNC_PROPERTIES_HH
#define DUMUX_2PNC_PROPERTIES_HH

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

//! The type tag for the implicit isothermal two phase n component problems
NEW_TYPE_TAG(TwoPNC);
NEW_TYPE_TAG(BoxTwoPNC, INHERITS_FROM(BoxModel, TwoPNC));
NEW_TYPE_TAG(CCTwoPNC, INHERITS_FROM(CCModel, TwoPNC));

//! The type tag for the implicit non-isothermal two phase n component problems
NEW_TYPE_TAG(TwoPNCNI, INHERITS_FROM(TwoPNC, NonIsothermal));
NEW_TYPE_TAG(BoxTwoPNCNI, INHERITS_FROM(BoxModel, TwoPNCNI));
NEW_TYPE_TAG(CCTwoPNCNI, INHERITS_FROM(CCModel, TwoPNCNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(NumMajorComponents); //!< Number of major fluid components which are considered in the calculation of the phase density
NEW_PROP_TAG(TwoPNCIndices); //!< Enumerations for the 2pncMin models
NEW_PROP_TAG(Formulation);   //!< The formulation of the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(Chemistry); //!< The chemistry class with which solves equlibrium reactions

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The parameters of the material law (extracted from the spatial parameters)

NEW_PROP_TAG(ReplaceCompEqIdx); //!< The index of the total mass balance equation, if one component balance is replaced (ReplaceCompEqIdx < NumComponents)
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the upwind weight for the mass conservation equations
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(BaseFluxVariables); //! The base flux variables
}
}

#endif
