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
 * \ingroup OnePNCModel
 *
 * \file
 *
 * \brief Defines the properties required for the one-phase n-component
 *        fully implicit model.
 */
#ifndef DUMUX_1PNC_PROPERTIES_HH
#define DUMUX_1PNC_PROPERTIES_HH

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

//! The type tag for the implicit isothermal one phase n component problems
NEW_TYPE_TAG(OnePNC);
NEW_TYPE_TAG(BoxOnePNC, INHERITS_FROM(BoxModel, OnePNC));
NEW_TYPE_TAG(CCOnePNC, INHERITS_FROM(CCModel, OnePNC));

//! The type tag for the implicit non-isothermal one phase n component problems
NEW_TYPE_TAG(OnePNCNI, INHERITS_FROM(OnePNC, NonIsothermal));
NEW_TYPE_TAG(BoxOnePNCNI, INHERITS_FROM(BoxModel, OnePNCNI));
NEW_TYPE_TAG(CCOnePNCNI, INHERITS_FROM(CCModel, OnePNCNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(PhaseIdx); //!< A phase index in to allow that a two-phase fluidsystem is used
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(FluidState); //!< Type of the fluid state to be used
//NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(UseMoles); //!< Defines whether mole (true) or mass (false) fractions are used
//NEW_PROP_TAG(BaseFluxVariables); //! The base flux variables
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(TauTortuosity); //!< Tortuosity value (tau) used in macroscopic diffusion

}
}

#endif
