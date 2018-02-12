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
 * \ingroup OnePTwoCAdsorptionModel
 * \file
 *
 * \brief Defines the properties required for the single-phase,
 *        two-component fully implicit model.
 */

#ifndef DUMUX_1P2C_ADSORPTION_PROPERTIES_HH
#define DUMUX_1P2C_ADSORPTION_PROPERTIES_HH


#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

namespace Dumux
{
// \{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
NEW_TYPE_TAG(OnePTwoCAdsorption);
NEW_TYPE_TAG(BoxOnePTwoCAdsorption, INHERITS_FROM(BoxModel, OnePTwoCAdsorption));
NEW_TYPE_TAG(CCOnePTwoCAdsorption, INHERITS_FROM(CCModel, OnePTwoCAdsorption));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(OnePTwoCAdsorptionNI, INHERITS_FROM(OnePTwoCAdsorption, NonIsothermal));
NEW_TYPE_TAG(BoxOnePTwoCAdsorptionNI, INHERITS_FROM(BoxModel, OnePTwoCAdsorptionNI));
NEW_TYPE_TAG(CCOnePTwoCAdsorptionNI, INHERITS_FROM(CCModel, OnePTwoCAdsorptionNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(PhaseIdx); //!< A phase index in to allow that a two-phase fluidsystem is used
NEW_PROP_TAG(NumComponents);   //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(FluidState); //!< Type of the fluid state to be used
NEW_PROP_TAG(ImplicitMassUpwindWeight);   //!< The default value of the upwind weight
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(UseMoles); //!< Defines whether mole (true) or mass (false) fractions are used
NEW_PROP_TAG(Scaling); //!< Defines Scaling of the model
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(TauTortuosity); //!< Tortuosity value (tau) used in macroscopic diffusion
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(SorptionCoefficientsuseAdsorption); //!< Determines if Adsorption in general is used
NEW_PROP_TAG(SorptionCoefficientsuseEL); //!< Determines if extended Langmuir Ad- and Desorption is used
NEW_PROP_TAG(SorptionCoefficientsuseL); //!< Determines if Langmuir Ad- and Desorption is used
NEW_PROP_TAG(SorptionCoefficientsuseF); //!< Determines if Freundlich Ad- and Desorption is used
NEW_PROP_TAG(SorptionCoefficientsuseBET); //!< Determines if Freundlich Ad- and Desorption is used


}
// \}
}

#endif
