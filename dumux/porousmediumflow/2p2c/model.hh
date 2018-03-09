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
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase two-component fully implicit model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{ D_{\alpha,\text{pm}}^\kappa \varrho_{\alpha} \frac{M^\kappa}{M_\alpha}
 \textbf{grad} x^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$x^\kappa_w + x^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPTwoCFormulation::pwsn</tt> or <tt>TwoPTwoCFormulation::pnsw</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 * Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^a_w<x^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^w_n<x^w_{n,max})\f$</li>
 * </ul>
 */
#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include <dune/common/fvector.hh>

// property forward declarations
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "primaryvariableswitch.hh"
#include "vtkoutputfields.hh"

namespace Dumux {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(TwoPTwoC, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(TwoPTwoCNI, INHERITS_FROM(TwoPTwoC, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the number of equations to 2
SET_INT_PROP(TwoPTwoC, NumEq, 2);

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use a static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPTwoC, NumComponents)
{
    static constexpr int value = 2;
    static_assert(GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents == value,
                  "Only fluid systems with 2 components are supported by the 2p-2c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use a static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPTwoC, NumPhases)
{
    static constexpr int value = 2;
    static_assert(GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases == value,
                  "Only fluid systems with 2 phases are supported by the 2p-2c model!");
};

//! Set the vtk output fields specific to this model
SET_PROP(TwoPTwoC, VtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = TwoPTwoCVtkOutputFields<FluidSystem, Indices>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(TwoPTwoC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the default formulation to pw-sn
SET_INT_PROP(TwoPTwoC, Formulation, TwoPTwoCFormulation::pwsn);

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(TwoPTwoC, ReplaceCompEqIdx, GET_PROP_VALUE(TypeTag, NumComponents));

//! Use the compositional local residual operator
SET_TYPE_PROP(TwoPTwoC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(TwoPTwoC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(TwoPTwoC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(TwoPTwoC, EnableEnergyBalance, false);

//! The primary variable switch for the 2p2c model
SET_TYPE_PROP(TwoPTwoC, PrimaryVariableSwitch, TwoPTwoCPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 2p2c model
SET_PROP(TwoPTwoC, PrimaryVariables)
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                     GET_PROP_VALUE(TypeTag, NumEq)>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Use the 2p2c VolumeVariables
SET_TYPE_PROP(TwoPTwoC, VolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//! Set the indices required by the isothermal 2p2c
SET_TYPE_PROP(TwoPTwoC, Indices, TwoPTwoCIndices<typename GET_PROP_TYPE(TypeTag, FluidSystem), /*PVOffset=*/0>);

//! Use the FVSpatialParams by default
SET_TYPE_PROP(TwoPTwoC, SpatialParams, FVSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(TwoPTwoC, EffectiveDiffusivityModel,
             DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(TwoPTwoC, UseMoles, true);

//! Determines whether the constraint solver is used
SET_BOOL_PROP(TwoPTwoC, UseConstraintSolver, true);

//! Determines whether the Kelvin equation is used to adapt the saturation vapor pressure
SET_BOOL_PROP(TwoPTwoC, UseKelvinEquation, false);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPTwoCNI, ThermalConductivityModel)
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

//! Set isothermal Indices
SET_TYPE_PROP(TwoPTwoCNI, IsothermalIndices, TwoPTwoCIndices<typename GET_PROP_TYPE(TypeTag, FluidSystem), /*PVOffset=*/0>);

//! Set isothermal output fields
SET_PROP(TwoPTwoCNI, IsothermalVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = TwoPTwoCVtkOutputFields<FluidSystem, Indices>;
};

// Set isothermal NumEq
SET_INT_PROP(TwoPTwoCNI, IsothermalNumEq, 2);

} // end namespace Properties
} // end namespace Dumux

#endif
