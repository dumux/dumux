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
 * \ingroup ThreePWaterOilModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the three-phase two-component
 *        fully implicit model. This model requires the use of the non-isothermal
 *        extension found in dumux/implicit/nonisothermal
 */
#ifndef DUMUX_3P2CNI_PROPERTIES_HH
#define DUMUX_3P2CNI_PROPERTIES_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(ThreePWaterOilNI, INHERITS_FROM(NonIsothermal));
NEW_TYPE_TAG(BoxThreePWaterOilNI, INHERITS_FROM(BoxModel, ThreePWaterOilNI));
NEW_TYPE_TAG(CCThreePWaterOilNI, INHERITS_FROM(CCModel, ThreePWaterOilNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)

NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem

NEW_PROP_TAG(UseMoles); //!Defines whether mole (true) or mass (false) fractions are used

NEW_PROP_TAG(UseMassOutput); //!Defines whether mole or mass are used for phaseStorage output
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity

NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(UseSimpleModel); //!< Determines whether a simple model with only two phase states (wng) and (wn) should be used
NEW_PROP_TAG(BaseFluxVariables); //! The base flux variables
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
}
}

#endif
