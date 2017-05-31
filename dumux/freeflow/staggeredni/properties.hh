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
#ifndef DUMUX_NAVIERSTOKES_NI_PROPERTIES_HH
#define DUMUX_NAVIERSTOKES_NI_PROPERTIES_HH

// #include <dumux/freeflow/staggered/properties.hh>

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

//! The type tags for the non-isothermal Navier Stokes problems
NEW_TYPE_TAG(StaggeredNonIsothermal);

NEW_PROP_TAG(EnableEnergyBalanceStokes);

//! The type tags for the corresponding non-isothermal problems
// NEW_TYPE_TAG(NavierStokesNI, INHERITS_FROM(NavierStokes, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(PhaseIdx); //!< Defines the phaseIdx


//////////////////////////////////////////////////////////////////
// Property tags required for the non-isothermal models
//////////////////////////////////////////////////////////////////

//TODO cleanup

NEW_PROP_TAG(IsothermalModel);
NEW_PROP_TAG(IsothermalFluxVariables);
NEW_PROP_TAG(IsothermalVolumeVariables);
NEW_PROP_TAG(IsothermalIndices);
NEW_PROP_TAG(IsothermalNumEqCellCenter);
NEW_PROP_TAG(IsothermalNumEqFace);
NEW_PROP_TAG(HaveVariableFormulation);
NEW_PROP_TAG(ThermalConductivityModel);
NEW_PROP_TAG(NiOutputLevel);

// forward declaration of other property tags
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(FluidSystem);


// \}
}

} // end namespace

#endif
