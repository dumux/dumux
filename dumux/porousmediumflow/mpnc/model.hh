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
 * \ingroup MPNCModel
 * \brief A fully implicit model for MpNc flow using
 *        vertex centered finite volumes.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 * The momentum approximation can be selected via "BaseFluxVariables":
 * Darcy (ImplicitDarcyFluxVariables) and Forchheimer (ImplicitForchheimerFluxVariables)
 * relations are available for all Box models.
 *
 * By inserting this into the equations for the conservation of the
 * mass of each component, one gets one mass-continuity equation for
 * each component \f$\kappa\f$
 * \f[
 \sum_{\kappa} \left(
    \phi \frac{\partial \left(\varrho_\alpha x_\alpha^\kappa S_\alpha\right)}{\partial t}
    +
    \mathrm{div}\;
    \left\{
        v_\alpha
       \frac{\varrho_\alpha}{\overline M_\alpha} x_\alpha^\kappa
    \right\}
    \right)
    = q^\kappa
    \f]
 * with \f$\overline M_\alpha\f$ being the average molar mass of the
 * phase \f$\alpha\f$: \f[ \overline M_\alpha = \sum_\kappa M^\kappa
 * \; x_\alpha^\kappa \f]
 *
 * For the missing \f$M\f$ model assumptions, the model assumes that
 * if a fluid phase is not present, the sum of the mole fractions of
 * this fluid phase is smaller than \f$1\f$, i.e.
 * \f[
 * \forall \alpha: S_\alpha = 0 \implies \sum_\kappa x_\alpha^\kappa \leq 1
 * \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be positive:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0 \f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[ \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0 \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$ which have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * These equations are then discretized using a fully-implicit vertex
 * centered finite volume scheme (often known as 'box'-scheme) for
 * spatial discretization and the implicit Euler method as temporal
 * discretization.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The pressure of the first phase \f$p_1\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */

#ifndef DUMUX_MPNC_MODEL_HH
#define DUMUX_MPNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysimplefluidlumping.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonequilibrium/model.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "localresidual.hh"

/*!
 * \ingroup \ingroup MPNCModel
 * \brief  Defines the properties required for the MpNc fully implicit model.
 */
namespace Dumux
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tags for the isothermal & non-isothermal two-phase model
NEW_TYPE_TAG(MPNC, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(MPNCNI, INHERITS_FROM(MPNC, NonIsothermal));
NEW_TYPE_TAG(MPNCNonequil, INHERITS_FROM(MPNC, NonEquilibrium));

/////////////////////////////////////////////////////////////////
// Properties for the isothermal mpnc model
//////////////////////////////////////////////////////////////////

//! The indices required by the mpnc model
SET_PROP(MPNC, Indices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numEquationBalance = GET_PROP_VALUE(TypeTag, NumEqBalance);
public:
    using type = MPNCIndices<FluidSystem, numEquationBalance, /*PVOffset=*/0>;
};
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(MPNC, SpatialParams, FVSpatialParams<TypeTag>);

//! Use the MpNc local residual for the MpNc model
SET_TYPE_PROP(MPNC,
              LocalResidual,
              MPNCLocalResidual<TypeTag>);
/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_PROP(MPNC, NumEq)
{
private:
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    static const unsigned int value = Indices::numPrimaryVars;
};
/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system.
 */
SET_PROP(MPNC, NumComponents)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const unsigned int value = FluidSystem::numComponents;
};

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_PROP(MPNC, NumPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const unsigned int value = FluidSystem::numPhases;
};

/*!
 * \brief Set the property for the number of equations without energyEquations.
 * In other models this would just be number of mass balance equations. but in Mpnc we have additional ncp-equations which also need to be added
 */
SET_INT_PROP(MPNC, NumEqBalance, GET_PROP_VALUE(TypeTag, NumPhases)+GET_PROP_VALUE(TypeTag, NumPhases));

//! the VolumeVariables property
SET_TYPE_PROP(MPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! This model uses the compositional fluid state
SET_PROP(MPNC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

SET_INT_PROP(MPNC, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents); //! Per default, no component mass balance is replaced

SET_BOOL_PROP(MPNC, UseMoles, true);                         //! Use mole fractions in the balance equations by default

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(MPNC, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);


//! Set the default pressure formulation to the pressure of the (most) wetting phase
SET_INT_PROP(MPNC,
             PressureFormulation,
             MpNcPressureFormulation::mostWettingFirst);

//! Set the vtk output fields specific to this model
SET_PROP(MPNC, VtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = MPNCVtkOutputFields<FluidSystem, Indices>;
};

SET_BOOL_PROP(MPNC, EnableAdvection, true);                  //! Enable advection
SET_BOOL_PROP(MPNC, EnableMolecularDiffusion, true);         //! Enable molecular diffusion
SET_BOOL_PROP(MPNC, EnableEnergyBalance, false);             //! This is the isothermal variant of the model

SET_BOOL_PROP(MPNC, EnableThermalNonEquilibrium, false);
SET_BOOL_PROP(MPNC, EnableChemicalNonEquilibrium, false);

SET_INT_PROP(MPNC, NumEnergyEqFluid, 0);
SET_INT_PROP(MPNC, NumEnergyEqSolid, 0);


/////////////////////////////////////////////////
// Properties for the non-isothermal mpnc model
/////////////////////////////////////////////////
SET_TYPE_PROP(MPNCNI, IsothermalVolumeVariables, MPNCVolumeVariables<TypeTag>);    //! set isothermal VolumeVariables
SET_TYPE_PROP(MPNCNI, IsothermalLocalResidual, MPNCLocalResidual<TypeTag>); //! set isothermal LocalResidual
//! set isothermal Indices
SET_PROP(MPNCNI, IsothermalIndices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numEquationBalance = GET_PROP_VALUE(TypeTag, NumEqBalance);
public:
    using type = MPNCIndices<FluidSystem, numEquationBalance, /*PVOffset=*/0>;
};
//! set isothermal NumEq
SET_INT_PROP(MPNCNI, IsothermalNumEq, GET_PROP_VALUE(TypeTag, NumEq));

/////////////////////////////////////////////////
// Properties for the non-equilibrium mpnc model
/////////////////////////////////////////////////

SET_TYPE_PROP(MPNCNonequil, EquilibriumLocalResidual, MPNCLocalResidual<TypeTag>);

//! Set the vtk output fields specific to this model
SET_PROP(MPNCNonequil, EquilibriumVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = MPNCVtkOutputFields<FluidSystem, Indices>;
};

SET_PROP(MPNCNonequil, EquilibriumIndices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numEquationBalance = GET_PROP_VALUE(TypeTag, NumEqBalance);
public:
    using type = MPNCIndices<FluidSystem, numEquationBalance, /*PVOffset=*/0>;
};

//number of balance equations means all transport equations and the constraint equations, energy equations come on top of that
SET_INT_PROP(MPNCNonequil, NumEqBalance, GET_PROP_VALUE(TypeTag, NumPhases)*GET_PROP_VALUE(TypeTag, NumComponents)+GET_PROP_VALUE(TypeTag, NumPhases));

//! in case we do not assume full non-equilibrium one needs a thermal conductivity
SET_PROP(MPNCNonequil, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySimpleFluidLumping<Scalar, GET_PROP_VALUE(TypeTag, NumEnergyEqFluid), Indices>;
};

} //end namespace Properties
} //end namespace Dumux

#endif
