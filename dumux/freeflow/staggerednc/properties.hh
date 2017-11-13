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
#ifndef DUMUX_NAVIERSTOKES_NC_PROPERTIES_HH
#define DUMUX_NAVIERSTOKES_NC_PROPERTIES_HH

#include <dumux/freeflow/staggered/properties.hh>
#include <dumux/freeflow/staggeredni/properties.hh>
#include <dumux/discretization/fickslaw.hh>

#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "vtkoutputfields.hh"

#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/material/fluidstates/compositional.hh>

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
NEW_TYPE_TAG(NavierStokesNC, INHERITS_FROM(NavierStokes));
NEW_TYPE_TAG(NavierStokesNCNI, INHERITS_FROM(NavierStokesNC, NavierStokesNonIsothermal));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
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
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    public:
        using type = CompositionalFluidState<Scalar, FluidSystem>;
};

SET_BOOL_PROP(NavierStokesNC, EnableComponentTransport, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(NavierStokesNC, EnableMolecularDiffusion, true);

SET_TYPE_PROP(NavierStokesNC, MolecularDiffusionType, FicksLaw<TypeTag>);

SET_BOOL_PROP(NavierStokesNC, UseMoles, false); //!< Defines whether molar (true) or mass (false) density is used

SET_INT_PROP(NavierStokesNC, PhaseIdx, 0); //!< Defines the phaseIdx

SET_TYPE_PROP(NavierStokesNC, VtkOutputFields, NavierStokesNCVtkOutputFields<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(NavierStokesNCNI, IsothermalIndices, NavierStokesNCIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(NavierStokesNCNI, IsothermalNumEqCellCenter, GET_PROP_VALUE(TypeTag, NumComponents)); //!< set the number of equations to 1
SET_INT_PROP(NavierStokesNCNI, IsothermalNumEqFace, 1); //!< set the number of equations


// \}
}

} // end namespace

#endif
