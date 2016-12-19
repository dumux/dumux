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
 * \ingroup NIModel
 * \file
 *
 * \brief Defines the properties required for the
 *        implicit non-isothermal models.
 */

#ifndef DUMUX_ENERGY_PROPERTIES_HH
#define DUMUX_ENERGY_PROPERTIES_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{
// \{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(NonIsothermal);

//////////////////////////////////////////////////////////////////
// Property tags required for the non-isothermal models
//////////////////////////////////////////////////////////////////

//TODO cleanup

NEW_PROP_TAG(IsothermalModel);
NEW_PROP_TAG(IsothermalFluxVariables);
NEW_PROP_TAG(IsothermalVolumeVariables);
NEW_PROP_TAG(IsothermalLocalResidual);
NEW_PROP_TAG(IsothermalIndices);
NEW_PROP_TAG(IsothermalNumEq);
NEW_PROP_TAG(HaveVariableFormulation);
NEW_PROP_TAG(ThermalConductivityModel);
NEW_PROP_TAG(NiOutputLevel);

// forward declaration of other property tags
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(FluidSystem);

}
// \}
}

#endif
