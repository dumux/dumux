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
 * \ingroup BoxStokesModel
 *
 * \file
 *
 * \brief Defines the properties required for the Stokes box model.
 */

#ifndef DUMUX_STOKESPROPERTIES_HH
#define DUMUX_STOKESPROPERTIES_HH

#include <dumux/implicit/box/properties.hh>


//added to resemble geomechanics
#include <dumux/implicit/fem/properties.hh>


namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the stokes problems
// deviated so it inherits from FemModel instead of BoxModel
NEW_TYPE_TAG(BoxStokes, INHERITS_FROM(FemModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(Fluid);
NEW_PROP_TAG(FluidSystem); //!< The employed fluid system
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(StokesStabilizationAlpha); //!< The parameter for the stabilization
NEW_PROP_TAG(StokesStabilizationBeta); //!< The parameter for the stabilization at boundaries
NEW_PROP_TAG(EnableUnsymmetrizedVelocityGradient); //!< Returns whether unsymmetrized velocity gradient for viscous term is used
NEW_PROP_TAG(EnableNavierStokes); //!< Returns whether Navier-Stokes should be solved instead of plain Stokes
NEW_PROP_TAG(EnablePseudoThreeDWallFriction); //!< Returns whether an additional wall friction term should be considered to mimic 3D behavior
NEW_PROP_TAG(UseMoles); //!< Defines whether molar (true) or mass (false) density is used

NEW_PROP_TAG(PhaseIdx); //!< A phase index in case that a two-phase fluidsystem is used
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(Scaling); //!<Defines Scaling of the model
}
}

#endif
