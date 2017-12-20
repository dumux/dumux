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
 * \file
 *
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_STAGGERED_NI_MODEL_HH
#define DUMUX_STAGGERED_NI_MODEL_HH

#include <dumux/common/properties.hh>
#include "indices.hh"
#include "vtkoutputfields.hh"
#include <dumux/discretization/fourierslaw.hh>


namespace Dumux
{

namespace Properties {

//! The type tags for the non-isothermal Navier Stokes problems
NEW_TYPE_TAG(NavierStokesNonIsothermal);

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////


SET_PROP(NavierStokesNonIsothermal, NumEq)
{
private:
    static constexpr auto isothermalNumEq = GET_PROP_VALUE(TypeTag, IsothermalNumEq);
public:
    static constexpr int value = isothermalNumEq + 1;
};

SET_TYPE_PROP(NavierStokesNonIsothermal, Indices, NavierStokesNonIsothermalIndices<TypeTag>);


SET_BOOL_PROP(NavierStokesNonIsothermal, EnableEnergyBalance, true);
SET_TYPE_PROP(NavierStokesNonIsothermal, VtkOutputFields, NavierStokesNonIsothermalVtkOutputFields<TypeTag>);

SET_TYPE_PROP(NavierStokesNonIsothermal, HeatConductionType, FouriersLaw<TypeTag>);

} // end namespace Properties

}

#endif
