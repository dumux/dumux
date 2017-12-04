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

#include <dumux/freeflow/properties.hh>

#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/freeflow/staggeredni/properties.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "indices.hh"
#include "velocityoutput.hh"
#include "vtkoutputfields.hh"
#include "boundarytypes.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/common/properties.hh>

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
NEW_TYPE_TAG(NavierStokes, INHERITS_FROM(FreeFlow));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(NavierStokesNI, INHERITS_FROM(NavierStokes, NavierStokesNonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(EnableInertiaTerms); //!< Returns whether to include inertia terms in the momentum balance eq or not (Stokes / Navier-Stokes)
NEW_PROP_TAG(EnableComponentTransport); //!< Returns whether to consider component transport or not
NEW_PROP_TAG(EnableEnergyTransport); //!<  Returns whether to consider energy transport or not
NEW_PROP_TAG(NormalizePressure); //!<  Returns whether to normalize the pressure term in the momentum balance or not
NEW_PROP_TAG(EnergyLocalResidual); //!<  The energy local residual
NEW_PROP_TAG(EnergyFluxVariables); //!<  The energy flux variables

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(NavierStokes, NumEqCellCenter, 1); //! set the number of equations to 1
SET_INT_PROP(NavierStokes, NumPhases, 1); //! The number of phases in the 1p model is 1
SET_INT_PROP(NavierStokes, NumComponents, 1); //! The number of components in the 1p model is 1

//! The fluid system to use by default
SET_TYPE_PROP(NavierStokes, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(NavierStokes, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokes, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> type;
};

//! The local residual function
SET_TYPE_PROP(NavierStokes, LocalResidual, StaggeredNavierStokesResidual<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(NavierStokes, VolumeVariables, NavierStokesVolumeVariables<TypeTag>);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(NavierStokes, FluxVariables, FreeFlowFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(NavierStokes, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

//! Enable advection
SET_BOOL_PROP(NavierStokes, EnableAdvection, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(NavierStokes, EnableMolecularDiffusion, false);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(NavierStokes, Indices, NavierStokesCommonIndices<TypeTag>);

SET_TYPE_PROP(NavierStokes, EnergyLocalResidual, FreeFlowEnergyLocalResidual<TypeTag>);

SET_TYPE_PROP(NavierStokes, EnergyFluxVariables, FreeFlowEnergyFluxVariables<TypeTag>);

SET_BOOL_PROP(NavierStokes, EnableEnergyBalance, false);

SET_TYPE_PROP(NavierStokes, VelocityOutput, StaggeredFreeFlowVelocityOutput<TypeTag>);

SET_TYPE_PROP(NavierStokes, VtkOutputFields, NavierStokesVtkOutputFields<TypeTag>);

SET_BOOL_PROP(NavierStokes, EnableInertiaTerms, true);

SET_BOOL_PROP(NavierStokes, EnableEnergyTransport, false);

SET_BOOL_PROP(NavierStokes, EnableComponentTransport, false);

//! Normalize the pressure term in the momentum balance or not
SET_BOOL_PROP(NavierStokes, NormalizePressure, true);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//set isothermal Indices
SET_TYPE_PROP(NavierStokesNI, IsothermalIndices, NavierStokesCommonIndices<TypeTag>);
SET_TYPE_PROP(NavierStokesNI, IsothermalVtkOutputFields, NavierStokesVtkOutputFields<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(NavierStokesNI, IsothermalNumEqCellCenter, 1); //!< set the number of equations to 1
SET_INT_PROP(NavierStokesNI, IsothermalNumEqFace, 1); //!< set the number of equations to 1

// \}
}

} // end namespace

#endif
