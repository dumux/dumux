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
 *  \file
 * \ingroup TwoPNCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 * v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 * \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}{
 * && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 * {\partial t}
 * - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 * \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 * (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mathbf{g}) \right\}
 * \nonumber \\ \nonumber \\
 *   &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 * - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 * \alpha \in \{w, g\}
 * \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as
 * spatial and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pwsn or TwoPTwoCIndices::pnsw. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 *Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only nonwetting phase is present: The mole fraction of, e.g., water in the nonwetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 */

#ifndef DUMUX_2PNC_MODEL_HH
#define DUMUX_2PNC_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

#include "volumevariables.hh"
#include "iofields.hh"
#include "indices.hh"

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief Specifies a number properties of two-phase n-component models.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use molar or mass balances
 * \tparam setMoleFractionForFP whether to set mole fractions for first or second phase
 */
template<int nComp, bool useMol, bool setMoleFractionForFP, TwoPFormulation formulation, int repCompEqIdx = nComp>
struct TwoPNCModelTraits
{
    using Indices = TwoPNCIndices;

    static constexpr int numEq() { return nComp; }
    static constexpr int numFluidPhases() { return 2; }
    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
    static constexpr bool enableThermalNonEquilibrium() { return false; }
    static constexpr bool enableChemicalNonEquilibrium() { return false; }

    static constexpr bool useMoles() { return useMol; }
    static constexpr bool setMoleFractionsForFirstPhase() { return setMoleFractionForFP; }

    static constexpr TwoPFormulation priVarFormulation() { return formulation; }
};

namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct TwoPNC { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct TwoPNCNI { using InheritsFrom = std::tuple<TwoPNC>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Properties for the isothermal 2pnc model
//////////////////////////////////////////////////////////////////
//! The primary variables vector for the 2pnc model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::TwoPNC>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

    template<class BaseTraits, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };

public:
    using type = TwoPNCVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

//! Set the base model traits
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::TwoPNC>
{
private:
    //! we use the number of components specified by the fluid system here
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numPhases == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
public:
    using type = TwoPNCModelTraits<FluidSystem::numComponents,
                                   getPropValue<TypeTag, Properties::UseMoles>(),
                                   getPropValue<TypeTag, Properties::SetMoleFractionsForFirstPhase>(),
                                   getPropValue<TypeTag, Properties::Formulation>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPNC> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPNC> { using type = TwoPNCIOFields; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPNC> { using type = CompositionalLocalResidual<TypeTag>; };                  //!< Use the compositional local residual

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::TwoPNC> { static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::numComponents; }; //!< Per default, no component mass balance is replaced

//! Default formulation is pw-Sn, overwrite if necessary
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPNC>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct SetMoleFractionsForFirstPhase<TypeTag, TTag::TwoPNC> { static constexpr bool value = true; };  //!< Set the primary variables mole fractions for the wetting or nonwetting phase
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPNC> { static constexpr bool value = true; };                         //!< Use mole fractions in the balance equations by default

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::TwoPNC> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! This model uses the compositional fluid state
template<class TypeTag>
struct FluidState<TypeTag, TTag::TwoPNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

/////////////////////////////////////////////////
// Properties for the non-isothermal 2pnc model
/////////////////////////////////////////////////

//! Set the non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPNCNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPNCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class DT, class EDM, class ETCM>
    struct NCNITraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
        using EffectiveThermalConductivityModel = ETCM;
    };
public:
    using type = TwoPNCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPNCNI> { using type = EnergyIOFields<TwoPNCIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPNCNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
