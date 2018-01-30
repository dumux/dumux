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
 * \ingroup NavierStokesModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_NAVIERSTOKES_PROPERTIES_HH
#define DUMUX_NAVIERSTOKES_PROPERTIES_HH

#include <dumux/freeflow/staggeredni/properties.hh>


namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit single-phase problems
NEW_TYPE_TAG(NavierStokes);

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(NavierStokesNI, INHERITS_FROM(NavierStokes, NavierStokesNonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(FluidSystem); //!< The type of the fluid system to use
NEW_PROP_TAG(Fluid); //!< The fluid used for the default fluid system
NEW_PROP_TAG(FluidState); //!< The type of the fluid state to use
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< Returns weight of the upwind cell when calculating fluxes
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(EnableInertiaTerms); //!< Returns whether to include inertia terms in the momentum balance eq or not (Stokes / Navier-Stokes)
NEW_PROP_TAG(BoundaryValues); //!< Type to set values on the boundary
NEW_PROP_TAG(EnableComponentTransport); //!< Returns whether to consider component transport or not
NEW_PROP_TAG(EnableEnergyTransport); //!<  Returns whether to consider energy transport or not
NEW_PROP_TAG(FaceVariables); //!<  Returns whether to consider energy transport or not
NEW_PROP_TAG(NormalizePressure); //!<  Returns whether to normalize the pressure term in the momentum balance or not
NEW_PROP_TAG(EnergyLocalResidual); //!<  The energy local residual
NEW_PROP_TAG(EnergyFluxVariables); //!<  The energy flux variables
NEW_PROP_TAG(PhaseIdx);
// \}
}

} // end namespace

#endif
