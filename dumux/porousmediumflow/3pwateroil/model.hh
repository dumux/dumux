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
 * \ingroup ThreePWaterOilModel
 * \brief Adaption of the fully implicit scheme to the three-phase three-component flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 * This model implements three-phase two-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to two components
 * \f$\kappa \in \{ water, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_\alpha x_\alpha^\kappa \mbox{\bf K}
 (\textbf{grad}\, p_\alpha - \varrho_\alpha \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_\text{pm}^\kappa \varrho_\alpha \frac{M^\kappa}{M_\alpha}
 \textbf{grad} x^\kappa_{\alpha} \right\}
 - q^\kappa = 0 \qquad \forall \kappa , \; \forall \alpha
 \f}
 *
 * Note that these balance equations are molar.
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations and
 * \f$x^w_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. Different cases can be distinguished:
 * <ul>
 *  <li> ... to be completed ... </li>
 * </ul>
 */

#ifndef DUMUX_3P2CNI_MODEL_HH
#define DUMUX_3P2CNI_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include "indices.hh"
#include "model.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "primaryvariableswitch.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePWaterOilModel
 * \brief Adaption of the fully implicit scheme to the three-phase two-component flow model.
 *
 */
namespace Properties {

NEW_TYPE_TAG(ThreePWaterOilNI, INHERITS_FROM(PorousMediumFlow, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(ThreePWaterOilNI, NumComponents)
{
 private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

 public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 3p2cni model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(ThreePWaterOilNI, NumPhases)
{
 private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p2cni model!");
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(ThreePWaterOilNI, FluidState){
    private:
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    public:
        using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreePWaterOilNI, LocalResidual, ThreePWaterOilLocalResidual<TypeTag>);

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(ThreePWaterOilNI, ReplaceCompEqIdx, GET_PROP_VALUE(TypeTag, NumComponents));

//! Enable advection
SET_BOOL_PROP(ThreePWaterOilNI, EnableAdvection, true);

//! disable molecular diffusion for the 3p model
SET_BOOL_PROP(ThreePWaterOilNI, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(ThreePWaterOilNI, EnableEnergyBalance, true);

//! The primary variable switch for the 3p3c model
SET_TYPE_PROP(ThreePWaterOilNI, PrimaryVariableSwitch, ThreePWaterOilPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 3p3c model
SET_PROP(ThreePWaterOilNI, PrimaryVariables)
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                     GET_PROP_VALUE(TypeTag, NumEq)>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Determines whether a constraint solver should be used explicitly
SET_BOOL_PROP(ThreePWaterOilNI, OnlyGasPhaseCanDisappear, true);

//! the VolumeVariables property
SET_TYPE_PROP(ThreePWaterOilNI, VolumeVariables, ThreePWaterOilVolumeVariables<TypeTag>);

//! The spatial parameters to be employed.
//! Use FVSpatialParams by default.
SET_TYPE_PROP(ThreePWaterOilNI, SpatialParams, FVSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(ThreePWaterOilNI, EffectiveDiffusivityModel,
             DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// disable velocity output by default

SET_BOOL_PROP(ThreePWaterOilNI, UseMoles, true); //!< Define that mole fractions are used in the balance equations per default

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePWaterOilNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySomerton<Scalar, Indices>;
};


//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! Set the isothermal vktoutputfields
SET_PROP(ThreePWaterOilNI, IsothermalVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = ThreePWaterOilVtkOutputFields<FluidSystem, Indices>;
};

//set isothermal Indices
SET_PROP(ThreePWaterOilNI, IsothermalIndices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ThreePWaterOilIndices<FluidSystem, /*PVOffset=*/0>;
};

//set isothermal NumEq
SET_INT_PROP(ThreePWaterOilNI, IsothermalNumEq, 2);

}

 }

#endif
