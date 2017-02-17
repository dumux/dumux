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
 * \ingroup OnePModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_NAVIER_STOKES_NC_PROPERTY_DEFAULTS_HH
#define DUMUX_NAVIER_STOKES_NC_PROPERTY_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "../staggered/problem.hh"
// #include "../staggered/model.hh"
#include "../staggered/propertydefaults.hh"


#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/material/fluidstates/compositional.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(FluxVariables);
NEW_PROP_TAG(FluxVariablesCache);
}
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

SET_PROP(NavierStokesNC, NumEqCellCenter)
{
private:
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
public:
    static constexpr int value = numComponents;
};


SET_INT_PROP(NavierStokesNC, ReplaceCompEqIdx, 0);


 /*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(NavierStokesNC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static constexpr int value = FluidSystem::numComponents;

};


//! the VolumeVariables property
SET_TYPE_PROP(NavierStokesNC, VolumeVariables, NavierStokesNCVolumeVariables<TypeTag>);
SET_TYPE_PROP(NavierStokesNC, Model, NavierStokesNCModel<TypeTag>);
SET_TYPE_PROP(NavierStokesNC, Indices, NavierStokesNCIndices<TypeTag>);


/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokesNC, FluidState)
{
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef CompositionalFluidState<Scalar, FluidSystem> type;
};

//
SET_BOOL_PROP(NavierStokesNC, EnableComponentTransport, true);

SET_BOOL_PROP(NavierStokesNC, UseMoles, false); //!< Defines whether molar (true) or mass (false) density is used

SET_INT_PROP(NavierStokesNC, PhaseIdx, 0); //!< Defines the phaseIdx


//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
// SET_TYPE_PROP(NavierStokesNI, IsothermalModel, NavierStokesModel<TypeTag>);

//set isothermal VolumeVariables
// SET_TYPE_PROP(NavierStokesNI, IsothermalVolumeVariables, NavierStokesVolumeVariables<TypeTag>);

//set isothermal LocalResidual
// SET_TYPE_PROP(NavierStokesNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
// SET_TYPE_PROP(NavierStokesNI, IsothermalIndices, NavierStokesCommonIndices<TypeTag>);

//set isothermal NumEq
// SET_INT_PROP(NavierStokesNI, IsothermalNumEq, 1);


// \}
} // end namespace Properties

} // end namespace Dumux

#endif
