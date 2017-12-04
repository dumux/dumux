// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup TwoPNCMinModel
 *
 * \file
 *
 * \brief Defines the properties required for the two-phase n-component mineralization
 *        fully implicit model.
 */
#ifndef DUMUX_2PNCMIN_PROPERTIES_HH
#define DUMUX_2PNCMIN_PROPERTIES_HH

#include <dumux/porousmediumflow/2pnc/implicit/properties.hh>

#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "localresidual.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(TwoPNCMin, INHERITS_FROM(TwoPNC));
NEW_TYPE_TAG(TwoPNCMinNI, INHERITS_FROM(TwoPNCMin, TwoPNCNI, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 2pncmin model
//////////////////////////////////////////////////////////////////

SET_TYPE_PROP(TwoPNCMin, VolumeVariables, TwoPNCMinVolumeVariables<TypeTag>);                  //! the VolumeVariables property
SET_TYPE_PROP(TwoPNCMin, VtkOutputFields, TwoPNCMinVtkOutputFields<TypeTag>);                  //! Set the vtk output fields specific to the TwoPNCMin model
SET_TYPE_PROP(TwoPNCMin, LocalResidual, TwoPNCMinLocalResidual<TypeTag>);                  //! Use the compositional local residual

//! Set the property for the number of solid phases, excluding the non-reactive matrix.
SET_PROP(TwoPNCMin, NumSPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));

public:
    static const int value = FluidSystem::numSPhases;
};

//! Set the property for the number of equations. For each component and each
//precipitated mineral/solid phase one equation has to be solved.

SET_PROP(TwoPNCMin, NumEq)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));

public:
    static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
};
/////////////////////////////////////////////////
// Properties for the non-isothermal 2pncmin model
/////////////////////////////////////////////////
SET_TYPE_PROP(TwoPNCMinNI, IsothermalVolumeVariables, TwoPNCMinVolumeVariables<TypeTag>);     //! set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNCMinNI, IsothermalVtkOutputFields, TwoPNCMinVtkOutputFields<TypeTag>);     //! set isothermal output fields

}
}

#endif
