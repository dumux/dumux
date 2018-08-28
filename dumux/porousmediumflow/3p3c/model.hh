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
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/vtkoutputfields.hh>

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
 * \param useCS if we are using the contraint solver
 */
template<bool useCS, bool useMol>
struct ThreePThreeCModelTraits
{
    using Indices = ThreePThreeCIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 3; }
    static constexpr int numComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool useConstraintSolver() { return useCS; }
    static constexpr bool useMoles() { return useMol; }

    template <class FluidSystem>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        switch (state)
        {
            case Indices::threePhases:
            {
                const std::vector<std::string> s1 = {"p_g",
                                                     "S_w",
                                                     "S_n"};
                return s1[pvIdx];
            }
            case Indices::wPhaseOnly:
            {
                const std::vector<std::string> s2 = {"p_g",
                                                     "x^" + FluidSystem::componentName(FluidSystem::gCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::wPhaseIdx),
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::wPhaseIdx)};
                return s2[pvIdx];
            }
            case Indices::gnPhaseOnly:
            {
                const std::vector<std::string> s3 = {"p_g",
                                                     "x^" + FluidSystem::componentName(FluidSystem::wCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx),
                                                     "S_n"};
                return s3[pvIdx];
            }
            case Indices::wnPhaseOnly:
            {
                const std::vector<std::string> s4 = {"p_g",
                                                     "x^" + FluidSystem::componentName(FluidSystem::gCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::wPhaseIdx),
                                                     "S_n"};
                return s4[pvIdx];
            }
            case Indices::gPhaseOnly:
            {
                const std::vector<std::string> s5 = {"p_g",
                                                     "x^" + FluidSystem::componentName(FluidSystem::wCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx),
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx)};
                return s5[pvIdx];
            }
            case Indices::wgPhaseOnly:
            {
                const std::vector<std::string> s6 = {"p_g",
                                                     "S_w",
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx)};
                return s6[pvIdx];
            }
        }
    }
};

/*!
 * \ingroup ThreePThreeCModel
 * \brief Traits class for the 3p3c model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct ThreePThreeCVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

namespace Properties {
//! The type tags for the isothermal three-phase three-component model
NEW_TYPE_TAG(ThreePThreeC, INHERITS_FROM(PorousMediumFlow));
//! The type tags for the non-isothermal three-phase three-component model
NEW_TYPE_TAG(ThreePThreeCNI, INHERITS_FROM(ThreePThreeC));

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
    using type = ThreePThreeCModelTraits<GET_PROP_VALUE(TypeTag, UseConstraintSolver), GET_PROP_VALUE(TypeTag, UseMoles)>;
};

//! Determines whether a constraint solver should be used explicitly
SET_BOOL_PROP(ThreePThreeC, UseConstraintSolver, false);

//! Set as default that _no_ component mass balance is replaced by the total mass balance
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
SET_TYPE_PROP(ThreePThreeC, PrimaryVariableSwitch, ThreePThreeCPrimaryVariableSwitch);

//! The primary variables vector for the 3p3c model
SET_PROP(ThreePThreeC, PrimaryVariables)
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                     GET_PROP_TYPE(TypeTag, ModelTraits)::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Set the volume variables property
SET_PROP(ThreePThreeC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    using Traits = ThreePThreeCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = ThreePThreeCVolumeVariables<Traits>;
};

//! The model after Millington (1961) is used for the effective diffusivity
SET_TYPE_PROP(ThreePThreeC, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(ThreePThreeC, VtkOutputFields, ThreePThreeCVtkOutputFields);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(ThreePThreeC, UseMoles, true);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_TYPE_PROP(ThreePThreeCNI, ThermalConductivityModel, ThermalConductivitySomerton<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! Set non-isothermal NumEq
SET_PROP(ThreePThreeCNI, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p3c model!");
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p3c model!");
    using IsothermalModelTraits = ThreePThreeCModelTraits<GET_PROP_VALUE(TypeTag, UseConstraintSolver), GET_PROP_VALUE(TypeTag, UseMoles)>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalModelTraits>;
};

//! Set the non-isothermal vktoutputfields
SET_TYPE_PROP(ThreePThreeCNI, VtkOutputFields, EnergyVtkOutputFields<ThreePThreeCVtkOutputFields>);

} // end namespace Properties
} // end namespace Dumux

#endif
