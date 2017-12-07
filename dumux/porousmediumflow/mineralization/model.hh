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
 * \ingroup Mineralization
 * \file
 *
 * \brief Defines the properties required for the
 *        implicit mineralization models.
 */

#ifndef DUMUX_MINERALIZATION_MODEL_HH
#define DUMUX_MINERALIZATION_MODEL_HH

#include <dumux/common/properties.hh>
#include "localresidual.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(Mineralization);

//! Set the general mineralization volume variables
SET_TYPE_PROP(Mineralization, VolumeVariables, MineralizationVolumeVariables<TypeTag>);

//! Set the general mineralization compositional local residual
SET_TYPE_PROP(Mineralization, LocalResidual, MineralizationLocalResidual<TypeTag>);

//! VTK outputs for mineralization models
SET_TYPE_PROP(Mineralization, VtkOutputFields, MineralizationVtkOutputFields<TypeTag>);

//! Set the property for the number of solid phases, excluding the non-reactive matrix.
SET_PROP(Mineralization, NumSPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static const int value = FluidSystem::numSPhases;
};

//! Set the property for the number of equations. For each component and each
//precipitated mineral/solid phase one equation has to be solved.
SET_PROP(Mineralization, NumEq)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
};

} // end namespace Properties
} // end namespace Dumux

#endif
