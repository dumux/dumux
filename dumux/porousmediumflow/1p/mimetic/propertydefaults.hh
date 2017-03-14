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
 * \ingroup OnepMimeticModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_MIMETIC_PROPERTY_DEFAULTS_HH
#define DUMUX_1P_MIMETIC_PROPERTY_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "../implicit/volumevariables.hh"
#include "indices.hh"
#include <dumux/porousmediumflow/immiscible/mimetic/localresidual.hh>
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "problem.hh"
#include "velocityoutput.hh"

#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(FluxVariables);
NEW_PROP_TAG(FluxVariablesCache);
NEW_PROP_TAG(StaggeredGeometryHelper);
}
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(OnePMimetic, NumEqCellCenter, 1); //!< set the number of equations to 1
SET_INT_PROP(OnePMimetic, NumEqFace, 1); //!< set the number of equations to 1
SET_INT_PROP(OnePMimetic, NumPhases, 1); //!< The number of phases in the 1p model is 1

//! The local residual function
SET_TYPE_PROP(OnePMimetic, LocalResidual, ImmiscibleMimeticLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(OnePMimetic, Model, OnePMimeticModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(OnePMimetic, VolumeVariables, OnePVolumeVariables<TypeTag>);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(OnePMimetic, FluxVariables, OnePMimeticFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(OnePMimetic, FluxVariablesCache, OnePMimeticFluxVariablesCache<TypeTag>);

//! Enable advection
SET_BOOL_PROP(OnePMimetic, EnableAdvection, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(OnePMimetic, EnableMolecularDiffusion, false);

//! Isothermal model by default
SET_BOOL_PROP(OnePMimetic, EnableEnergyBalance, false);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(OnePMimetic, Indices, OnePMimeticIndices<TypeTag>);


//! The fluid system to use by default
SET_TYPE_PROP(OnePMimetic, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(OnePMimetic, Fluid)
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
SET_PROP(OnePMimetic, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default
SET_BOOL_PROP(OnePMimetic, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(OnePMimetic, ProblemEnableGravity, true);

SET_BOOL_PROP(OnePMimetic, EnableInertiaTerms, false);

SET_BOOL_PROP(OnePMimetic, EnableEnergyTransport, false);

SET_BOOL_PROP(OnePMimetic, EnableComponentTransport, false);

SET_PROP(OnePMimetic, BoundaryValues)
{
private:
    using CellCenterBoundaryValues = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceBoundaryValues = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterBoundaryValues, FaceBoundaryValues>;
};

//SET_TYPE_PROP(NavierStokes, VtkOutputModule, FreeFlowStaggeredVtkOutputModule<TypeTag>);
//
SET_TYPE_PROP(OnePMimetic, VelocityOutput, OnePMimeticVelocityOutput<TypeTag>);

//! average is used as default model to compute the effective thermal heat conductivity
// SET_PROP(NavierStokesNI, ThermalConductivityModel)
// { private :
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   public:
//     typedef ThermalConductivityAverage<Scalar> type;
// };

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
