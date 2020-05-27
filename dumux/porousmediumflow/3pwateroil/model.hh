// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_\alpha x_\alpha^\kappa \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_\alpha \mathbf{g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_\text{pm}^\kappa \varrho_\alpha \frac{1}{M_\kappa}
 \textbf{grad} X^\kappa_{\alpha} \right\}
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
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include "indices.hh"
#include "model.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Specifies a number properties of the three-phase two-component model.
 */
template<bool onlyGasPhase>
struct ThreePWaterOilModelTraits
{
    using Indices = ThreePWaterOilIndices;

    static constexpr int numEq() { return 2; }
    static constexpr int numFluidPhases() { return 3; }
    static constexpr int numFluidComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool onlyGasPhaseCanDisappear() { return onlyGasPhase; }
};

namespace Properties {

// Create new type tags
namespace TTag {
struct ThreePWaterOilNI { using InheritsFrom = std::tuple<PorousMediumFlow>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the non-isothermal model traits property
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ThreePWaterOilNI>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numComponents == 2, "Only fluid systems with 2 components are supported by the 3p2cni model!");
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p2cni model!");
public:
    using type = PorousMediumFlowNIModelTraits<ThreePWaterOilModelTraits<getPropValue<TypeTag, Properties::OnlyGasPhaseCanDisappear>()>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::ThreePWaterOilNI>{
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    public:
        using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! The local residual function of the conservation equations
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ThreePWaterOilNI> { using type = ThreePWaterOilLocalResidual<TypeTag>; };

//! Set as default that no component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::ThreePWaterOilNI> { static constexpr int value = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents(); };

//! The primary variables vector for the 3p3c model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ThreePWaterOilNI>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Determines whether a constraint solver should be used explicitly
template<class TypeTag>
struct OnlyGasPhaseCanDisappear<TypeTag, TTag::ThreePWaterOilNI> { static constexpr bool value = true; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ThreePWaterOilNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using BaseTraits = ThreePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class DT, class EDM, class ETCM>
    struct NCNITraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
        using EffectiveThermalConductivityModel = ETCM;
    };

public:
    using type = ThreePWaterOilVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::ThreePWaterOilNI>
{ using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

// Define that mole fractions are used in the balance equations per default
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ThreePWaterOilNI> { static constexpr bool value = true; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::ThreePWaterOilNI> { using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the non-isothermal vkt output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::ThreePWaterOilNI> { using type = EnergyIOFields<ThreePWaterOilIOFields>; };

} // end namespace Properties
} // end namespace Dumux

#endif
