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
 * \ingroup TwoPNCMinModel, Mineralization
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model with additional solid/mineral phases.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 * Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 * The primary variable of the solid phases is the volume fraction \f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
 *
 * The source an sink terms link the mass balances of the n-transported component to the solid phases.
 * The porosity \f$\phi\f$ is updated according to the reduction of the initial (or solid-phase-free porous medium) porosity \f$\phi_0\f$
 * by the accumulated volume fractions of the solid phases:
 * \f$ \phi = \phi_0 - \sum (\phi_\lambda)\f$
 * Additionally, the permeability is updated depending on the current porosity.
 */
#ifndef DUMUX_2PNCMIN_MODEL_HH
#define DUMUX_2PNCMIN_MODEL_HH

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/2pnc/volumevariables.hh>

#include <dumux/material/solidstates/compositionalsolidstate.hh>

#include <dumux/porousmediumflow/mineralization/model.hh>
#include <dumux/porousmediumflow/mineralization/localresidual.hh>
#include <dumux/porousmediumflow/mineralization/volumevariables.hh>
#include <dumux/porousmediumflow/mineralization/iofields.hh>

#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct TwoPNCMin { using InheritsFrom = std::tuple<TwoPNC>; };
struct TwoPNCMinNI { using InheritsFrom = std::tuple<TwoPNCMin>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 2pncmin model
//////////////////////////////////////////////////////////////////

// use the mineralization local residual
SET_TYPE_PROP(TwoPNCMin, LocalResidual, MineralizationLocalResidual<TypeTag>);

//! use the mineralization volume variables together with the 2pnc vol vars
SET_PROP(TwoPNCMin, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    using Traits = TwoPNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
    using NonMinVolVars = TwoPNCVolumeVariables<Traits>;
public:
    using type = MineralizationVolumeVariables<Traits, NonMinVolVars>;
};

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(TwoPNCMin, IOFields, MineralizationIOFields<TwoPNCIOFields>);

//! The 2pnc model traits define the non-mineralization part
SET_PROP(TwoPNCMin, ModelTraits)
{
private:
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using NonMineralizationTraits = typename GET_PROP_TYPE(TypeTag, BaseModelTraits);
public:
    using type = MineralizationModelTraits<NonMineralizationTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;
};

//! The two-phase model uses the immiscible fluid state
SET_PROP(TwoPNCMin, SolidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
public:
    using type = CompositionalSolidState<Scalar, SolidSystem>;
};

//////////////////////////////////////////////////////////////////
// Properties for the non-isothermal 2pncmin model
//////////////////////////////////////////////////////////////////

//! Set non-isothermal model traits
SET_PROP(TwoPNCMinNI, ModelTraits)
{
private:
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using TwoPNCTraits = typename GET_PROP_TYPE(TypeTag, BaseModelTraits);
    using IsothermalTraits = MineralizationModelTraits<TwoPNCTraits, SolidSystem::numComponents, SolidSystem::numInertComponents>;
public:
    // the mineralization traits, based on 2pnc traits, are the isothermal traits
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! non-isothermal vtkoutput
SET_PROP(TwoPNCMinNI, IOFields)
{
    using MineralizationIOF = MineralizationIOFields<TwoPNCIOFields>;
    using type = EnergyIOFields<MineralizationIOF>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
