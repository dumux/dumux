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
 * \ingroup OnePNCNCModel
 *
 * \file
 *
 * \brief Defines the properties required for the one-phase n-component
 *        fully implicit model.
 */
#ifndef DUMUX_1PNC_PROPERTIES_HH
#define DUMUX_1PNC_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the implicit the isothermal & non-isothermal one phase n component problems
NEW_TYPE_TAG(OnePNC, INHERITS_FROM(PorousMediumFlow, NumericModel, LinearSolverTypeTag));
NEW_TYPE_TAG(OnePNCNI, INHERITS_FROM(OnePNC, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(OnePNC, NumComponents)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));

public:
    static constexpr auto value = FluidSystem::numComponents;
};

/*!
 * \brief Set the property for the number of equations: For each existing component one equation has to be solved.
 */
SET_PROP(OnePNC, NumEq)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static constexpr auto value = FluidSystem::numComponents;
};

//! Set as default that no component mass balance is replaced by the total mass balance
SET_PROP(OnePNC, ReplaceCompEqIdx)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static constexpr auto value = FluidSystem::numComponents;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OnePNC, FluidState){
    private:
        using Scalar =  typename GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
    public:
        using type = Dumux::CompositionalFluidState<Scalar, FluidSystem>;
};


SET_TYPE_PROP(OnePNC, LocalResidual, CompositionalLocalResidual<TypeTag>); //! The local residual function
SET_TYPE_PROP(OnePNC, VolumeVariables, OnePNCVolumeVariables<TypeTag>);   //! the VolumeVariables property
SET_BOOL_PROP(OnePNC, EnableAdvection, true);                           //! The one-phase model considers advection
SET_BOOL_PROP(OnePNC, EnableMolecularDiffusion, true);                 //! The one-phase model has no molecular diffusion
SET_BOOL_PROP(OnePNC, EnableEnergyBalance, false);                      //! Isothermal model by default
SET_TYPE_PROP(OnePNC, Indices, OnePNCIndices <TypeTag, /*PVOffset=*/0>);                            //! The indices required by the isothermal single-phase model
SET_TYPE_PROP(OnePNC, VtkOutputFields, OnePNCVtkOutputFields<TypeTag>);   //! Set the vtk output fields specific to this model


///////////////////////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of equations: For each existing component one equation has to be solved.
 */
SET_PROP(OnePNCNI, IsothermalNumEq)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static constexpr auto value = FluidSystem::numComponents;
};
SET_BOOL_PROP(OnePNCNI, EnableEnergyBalance, true);                                   //! we do solve for the energy balance here
SET_TYPE_PROP(OnePNCNI, IsothermalVtkOutputFields, OnePNCVtkOutputFields<TypeTag>);     //! the isothermal vtk output fields
SET_TYPE_PROP(OnePNCNI, IsothermalVolumeVariables, OnePNCVolumeVariables<TypeTag>);     //! Vol vars of the isothermal model
SET_TYPE_PROP(OnePNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);   //! Local residual of the isothermal model
SET_TYPE_PROP(OnePNCNI, IsothermalIndices, OnePNCIndices <TypeTag, /*PVOffset=*/0>);                              //! Indices of the isothermal model
SET_TYPE_PROP(OnePNCNI,
              ThermalConductivityModel,
              ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>); //! Use the average for effective conductivities

}
}

#endif
