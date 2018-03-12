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
 * \ingroup ThreePThreeCModel
 * \brief Adaption of the fully implicit scheme to the three-phase three-component
 *        flow model.
 *
 * This model implements three-phase three-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to three components
 * \f$\kappa \in \{ water, air, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_{\alpha,mol} x_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_{\alpha,mol} x_\alpha^\kappa \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_{\alpha,mass} \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_\text{pm}^\kappa \varrho_{\alpha,mol}
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
 * \f$x^w_\alpha + x^a_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases.
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. The following cases can be distinguished:
 * <ul>
 *  <li> All three phases are present: Primary variables are two saturations \f$(S_w\f$ and \f$S_n)\f$,
 *       and a pressure, in this case \f$p_g\f$. </li>
 *  <li> Only the water phase is present: Primary variables are now the mole fractions of air and
 *       contaminant in the water phase \f$(x_w^a\f$ and \f$x_w^c)\f$, as well as the gas pressure, which is,
 *       of course, in a case where only the water phase is present, just the same as the water pressure. </li>
 *  <li> Gas and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_g^w\f$, \f$p_g)\f$. </li>
 *  <li> Water and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_w^a\f$, \f$p_g)\f$. </li>
 *  <li> Only gas phase is present: Primary variables \f$(x_g^w\f$, \f$x_g^c\f$, \f$p_g)\f$. </li>
 *  <li> Water and gas phases are present: Primary variables \f$(S_w\f$, \f$x_w^g\f$, \f$p_g)\f$. </li>
 * </ul>
 */

#ifndef DUMUX_3P3C_MODEL_HH
#define DUMUX_3P3C_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "primaryvariableswitch.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief Specifies a number properties of two-phase models.
 */
struct ThreePThreeCModelTraits
{
    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 3; }
    static constexpr int numComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
};

namespace Properties {
//! The type tags for the isothermal three-phase three-component model
NEW_TYPE_TAG(ThreePThreeC, INHERITS_FROM(PorousMediumFlow));
//! The type tags for the non-isothermal three-phase three-component model
NEW_TYPE_TAG(ThreePThreeCNI, INHERITS_FROM(ThreePThreeC, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits
SET_PROP(ThreePThreeC, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p3c model!");
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p3c model!");
public:
    using type = ThreePThreeCModelTraits;
};

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(ThreePThreeC, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents());
/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(ThreePThreeC, FluidState){
    private:
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    public:
        using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreePThreeC, LocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//! The primary variable switch for the 3p3c model
SET_TYPE_PROP(ThreePThreeC, PrimaryVariableSwitch, ThreePThreeCPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 3p3c model
SET_PROP(ThreePThreeC, PrimaryVariables)
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                     GET_PROP_TYPE(TypeTag, ModelTraits)::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! the VolumeVariables property
SET_TYPE_PROP(ThreePThreeC, VolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//! Determines whether a constraint solver should be used explicitly
SET_BOOL_PROP(ThreePThreeC, UseConstraintSolver, false);

//! The indices required by the isothermal 3p3c model
SET_PROP(ThreePThreeC, Indices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ThreePThreeCIndices<FluidSystem, /*PVOffset=*/0>;
};

//! The spatial parameters to be employed.
//! Use FVSpatialParams by default.
SET_TYPE_PROP(ThreePThreeC, SpatialParams, FVSpatialParams<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(ThreePThreeC, EffectiveDiffusivityModel)
{ private :
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
 public:
    using type = DiffusivityMillingtonQuirk<Scalar>;
};

//! Set the vtk output fields specific to this model
SET_PROP(ThreePThreeC, VtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThreePThreeCVtkOutputFields<FluidSystem, Indices>;
};

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(ThreePThreeC, UseMoles, true);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePThreeCNI, ThermalConductivityModel)
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

//! Set isothermal VolumeVariables
SET_TYPE_PROP(ThreePThreeCNI, IsothermalVolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//! Set isothermal LocalResidual
SET_TYPE_PROP(ThreePThreeCNI, IsothermalLocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//! Set isothermal Indices
SET_PROP(ThreePThreeCNI, IsothermalIndices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ThreePThreeCIndices<FluidSystem, /*PVOffset=*/0>;
};

//! Set isothermal NumEq
SET_PROP(ThreePThreeCNI, IsothermalModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p3c model!");
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p3c model!");
public:
    using type = ThreePThreeCModelTraits;
};

//! Set the isothermal vktoutputfields
SET_PROP(ThreePThreeCNI, IsothermalVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = ThreePThreeCVtkOutputFields<FluidSystem, Indices>;
};

} // end namespace Properties

} // end namespace Dumux

#endif
