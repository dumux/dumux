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
 * \brief Adaption of the fully implicit scheme to the two-phase one-component flow model.
 *
 */
#ifndef DUMUX_2P1C_MODEL_HH
#define DUMUX_2P1C_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidstates/compositional.hh>


#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/model.hh>

#include <dumux/porousmediumflow/2p/implicit/vtkoutputfields.hh>

#include "darcyslaw.hh"
#include "vtkoutputfields.hh"
#include "localresidual.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "primaryvariableswitch.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPOneCModel
 * \brief Adaption of the fully implicit scheme to the two-phase one-component flow model.
 *
 * \note The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal.
 *
 * This model is designed for simulating two fluid phases with water as the only component.
 * It is particularly suitable for the simulation of steam injection in saturated conditions.
 *
 * The model implements the flow of two phases and one component, i.e. a pure liquid (e.g. water)
 * and its vapor (e.g. steam),
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
\phi \frac{\partial\ \sum_\alpha (\rho_\alpha S_\alpha)}{\partial t} \\-\sum \limits_ \alpha \text{div} \left \{\rho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}
\mathbf{K} (\mathbf{grad}p_\alpha - \rho_\alpha \mathbf{g}) \right \} -q^w =0
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. The model features a primary variable switch.
 * If only one phase is present, \f$p_g\f$ and \f$T\f$ are the primary variables.
 * In the presence of two phases, \f$p_g\f$ and \f$S_w\f$ become the new primary variables.
 */

namespace Properties
{
 //////////////////////////////////////////////////////////////////
 // Type tags
 //////////////////////////////////////////////////////////////////
 NEW_TYPE_TAG(TwoPOneCNI, INHERITS_FROM(PorousMediumFlow, NonIsothermal));

 //////////////////////////////////////////////////////////////////
 // Property tags
 //////////////////////////////////////////////////////////////////

 NEW_PROP_TAG(UseBlockingOfSpuriousFlow); //!< Determines whether Blocking ofspurious flow is used

 /*!
  * \brief Set the property for the number of components.
  *
  * We just forward the number from the fluid system and use an static
  * assert to make sure it is 1.
  */
 SET_PROP(TwoPOneCNI, NumComponents)
 {
  private:
     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
  public:
     static constexpr auto value = FluidSystem::numComponents;

     static_assert(value == 1,
                   "Only fluid systems with 1 component are supported by the 2p1cni model!");
 };

 /*!
  * \brief Set the property for the number of fluid phases.
  *
  * We just forward the number from the fluid system and use an static
  * assert to make sure it is 2.
  */
 SET_PROP(TwoPOneCNI, NumPhases)
 {
  private:
     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
  public:
     static constexpr auto value = FluidSystem::numPhases;
     static_assert(value == 2,
                   "Only fluid systems with 2 phases are supported by the 2p1cni model!");
 };

SET_TYPE_PROP(TwoPOneCNI, LocalResidual, TwoPOneCLocalResidual<TypeTag>); //! The local residual function

SET_BOOL_PROP(TwoPOneCNI, EnableAdvection, true);                           //! The one-phase model considers advection
SET_BOOL_PROP(TwoPOneCNI, EnableMolecularDiffusion, false);                 //! The one-phase model has no molecular diffusion

 /*!
  * \brief The fluid state which is used by the volume variables to
  *        store the thermodynamic state. This should be chosen
  *        appropriately for the model ((non-)isothermal, equilibrium, ...).
  *        This can be done in the problem.
  */
 SET_PROP(TwoPOneCNI, FluidState)
 {
 private:
     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
 public:
     using type = CompositionalFluidState<Scalar, FluidSystem>;
    //  using type = ImmiscibleFluidState<Scalar, FluidSystem>;
 };

 //! Determines whether Blocking ofspurious flow is used
 SET_BOOL_PROP(TwoPOneCNI, UseBlockingOfSpuriousFlow, false);

 SET_TYPE_PROP(TwoPOneCNI, AdvectionType, TwoPOneCDarcysLaw<TypeTag>);

  SET_TYPE_PROP(TwoPOneCNI, VolumeVariables, TwoPOneCVolumeVariables<TypeTag>);

 //! The primary variable switch for the 2p1c model
 SET_TYPE_PROP(TwoPOneCNI, PrimaryVariableSwitch, TwoPOneCPrimaryVariableSwitch<TypeTag>);

 //! The primary variables vector for the 2p1c model
 SET_TYPE_PROP(TwoPOneCNI, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);

 //! Somerton is used as default model to compute the effective thermal heat conductivity
 SET_PROP(TwoPOneCNI, ThermalConductivityModel)
 {
 private:
     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
     using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
 public:
     using type = ThermalConductivitySomerton<Scalar>;
 };

 //////////////////////////////////////////////////////////////////
 // Property values for isothermal model required for the general non-isothermal model
 //////////////////////////////////////////////////////////////////
 //set isothermal VolumeVariables
 SET_TYPE_PROP(TwoPOneCNI, IsothermalVolumeVariables, TwoPOneCVolumeVariables<TypeTag>);

 //set isothermal LocalResidual
 SET_TYPE_PROP(TwoPOneCNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

 //set isothermal Indices
 SET_TYPE_PROP(TwoPOneCNI, IsothermalIndices, TwoPOneCIndices<TypeTag, 0>);

//! the isothermal vtk output fields
 SET_TYPE_PROP(TwoPOneCNI, IsothermalVtkOutputFields, TwoPOneCVtkOutputFields<TypeTag>);

 //set isothermal NumEq
 SET_INT_PROP(TwoPOneCNI, IsothermalNumEq, 1);

} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_2P1C_MODEL_HH
