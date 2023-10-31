// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCModel
 * \brief Properties for a two-phase, two-component model for flow in porous media.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ \kappa_w, \kappa_n \}\f$, where \f$\kappa_w\f$ and \f$\kappa_n\f$ are
 * the main components of the wetting and nonwetting phases, respectively.
 * The governing equations are the mass or the mole conservation equations of the two components,
 * depending on the property <tt>UseMoles</tt>. The mass balance equations are given as
 * \f[
   \frac{\partial (\sum_\alpha \phi \rho_\alpha X_\alpha^\kappa S_\alpha)}{\partial t}
   - \sum_\alpha \nabla \cdot \left\{ \rho_\alpha X_\alpha^\kappa v_\alpha \right\}
   - \sum_\alpha \nabla \cdot \mathbf{F}_{\mathrm{diff, mass}, \alpha}^\kappa
   - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{\kappa_w, \kappa_n\} \, , \alpha \in \{w, n\}.
   \f]
 * The mole balance is given as
 * \f[
   \frac{\partial (\sum_\alpha \phi \varrho_{m, \alpha} x_\alpha^\kappa S_\alpha)}{\partial t}
   + \sum_\alpha \nabla \cdot \left\{ \varrho_{m, \alpha} x_\alpha^\kappa v_\alpha \right\}
   + \sum_\alpha \nabla \cdot \mathbf{F}_{\mathrm{diff, mole}, \alpha}^\kappa
   - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{\kappa_w, \kappa_n\} \, , \alpha \in \{w, n\},
   \f]
    * where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_\alpha \f$ represents the saturation of phase \f$ \alpha \f$,
 * * \f$ \rho_\alpha \f$ is the mass density of phase \f$ \alpha \f$,
 * * \f$ X_\alpha^\kappa \f$ is the mass fraction of component \f$ \kappa \f$ in phase  \f$ \alpha \f$,
 * * \f$ x_\alpha^\kappa \f$ is the mole fraction of component \f$ \kappa \f$ in phase    \f$ \alpha \f$,
 * * \f$ v_\alpha \f$ is the velocity of phase \f$ \alpha \f$,
 * * \f$ \mathbf{F}_{\mathrm{diff, mass}, \alpha}^\kappa \f$ represents the diffusive  mass flux of component \f$ \kappa \f$ in phase \f$ \alpha \f$,
 * * \f$ q_\alpha^\kappa \f$ is a source or sink term.
 *
 * Boundary conditions and sources have to be defined by the user in the corresponding
 * units. The default setting for the property <tt>UseMoles</tt> can be found in the 2pnc model.
 *
 * Per default, the Darcy's and Fick's law are used for the fluid phase velocities and the
 * diffusive fluxes, respectively. See dumux/flux/darcyslaw.hh and dumux/flux/fickslaw.hh
 * for more details.
 *
 * By using constitutive relations for the capillary pressure \f$p_c = p_n - p_w\f$ and
 * relative permeability \f$k_{r\alpha}\f$ and taking advantage of the fact that \f$S_w + S_n = 1\f$
 * and \f$x^{\kappa_w}_\alpha + x^{\kappa_n}_\alpha = 1\f$, the number of unknowns can be reduced to two.
 * In single-phase regimes, the used primary variables are either \f$p_w\f$ and \f$S_n\f$ (default)
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be specified by setting
 * the <tt>Formulation</tt> property to either
 * <tt>TwoPTwoCFormulation::pwsn</tt> or <tt>TwoPTwoCFormulation::pnsw</tt>.
 *
 * In two-phase flow regimes the second primary variable depends on the phase state and is the mole or mass
 * fraction (depending on the property <tt>UseMoles</tt>). The following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of the nonwetting phase main component in the wetting phase \f$x^{\kappa_n}_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^{\kappa_n}_w<x^{\kappa_n}_{w,max})\f$</li>
 *  <li> Only nonwetting phase is present: The mole fraction of the wetting phase main component in the nonwetting phase, \f$x^{\kappa_w}_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^{\kappa_w}_n<x^{\kappa_w}_{n,max})\f$</li>
 * </ul>
 */

#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include <array>

// property forward declarations
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/nonequilibrium/model.hh>
#include <dumux/porousmediumflow/nonequilibrium/volumevariables.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/simplefluidlumping.hh>

#include "volumevariables.hh"

namespace Dumux {

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
// Create new type tags
namespace TTag {
struct TwoPTwoC { using InheritsFrom = std::tuple<TwoPNC>; };
struct TwoPTwoCNI { using InheritsFrom = std::tuple<TwoPTwoC>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the model traits property.
 */
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::TwoPTwoC>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numComponents == 2, "Only fluid systems with 2 components are supported by the 2p-2c model!");
    static_assert(FluidSystem::numPhases == 2, "Only fluid systems with 2 phases are supported by the 2p-2c model!");

public:
    using type = TwoPNCModelTraits<FluidSystem::numComponents,
                                   getPropValue<TypeTag, Properties::UseMoles>(),
                                   /*setMFracForFirstPhase=*/true,
                                   getPropValue<TypeTag, Properties::Formulation>(),
                                   getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPTwoC> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; };

//! Use the 2p2c VolumeVariables
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    static_assert(FSY::numComponents == 2, "Only fluid systems with 2 components are supported by the 2p2c model!");
    static_assert(FSY::numPhases == 2, "Only fluid systems with 2 phases are supported by the 2p2c model!");
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
    static constexpr bool useConstraintSolver = getPropValue<TypeTag, Properties::UseConstraintSolver>();
public:
    using type = TwoPTwoCVolumeVariables<NCTraits<BaseTraits, DT, EDM>, useConstraintSolver>;
};

//! Determines whether the constraint solver is used
template<class TypeTag>
struct UseConstraintSolver<TypeTag, TTag::TwoPTwoC> { static constexpr bool value = true; };

//////////////////////////////////////////////////////////////////////
// Properties for the non-isothermal 2p2c model (inherited from 2pnc)
//////////////////////////////////////////////////////////////////////

//! Set the non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPTwoCNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
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
    using type = TwoPTwoCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>>;
};


//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPTwoCNI> { using type = EnergyIOFields<TwoPNCIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPTwoCNI> { using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Properties

template<class TwoPTwoCModelTraits>
struct TwoPTwoCUnconstrainedModelTraits : public TwoPTwoCModelTraits
{
    static constexpr int numConstraintEq() { return 0; }
};

namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
namespace TTag {
struct TwoPTwoCNonEquil { using InheritsFrom = std::tuple<NonEquilibrium, TwoPTwoC>; };
} // end namespace TTag

/////////////////////////////////////////////////
// Properties for the non-equilibrium TwoPTwoC model
/////////////////////////////////////////////////

template<class TypeTag>
struct EquilibriumLocalResidual<TypeTag, TTag::TwoPTwoCNonEquil> { using type = CompositionalLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct EquilibriumIOFields<TypeTag, TTag::TwoPTwoCNonEquil> { using type = TwoPNCIOFields; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPTwoCNonEquil>
{
private:
    using EquiTraits = GetPropType<TypeTag, Properties::EquilibriumModelTraits>;
    static constexpr bool enableTNE = getPropValue<TypeTag, Properties::EnableThermalNonEquilibrium>();
    static constexpr bool enableCNE = getPropValue<TypeTag, Properties::EnableChemicalNonEquilibrium>();
    static constexpr int numEF = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr int numES = getPropValue<TypeTag, Properties::NumEnergyEqSolid>();
    static constexpr auto nf = getPropValue<TypeTag, Properties::NusseltFormulation>();
    static constexpr auto ns = getPropValue<TypeTag, Properties::SherwoodFormulation>();

    using NonEquilTraits = NonEquilibriumModelTraits<EquiTraits, enableCNE, enableTNE, numEF, numES, nf, ns>;
public:
    using type = NonEquilTraits;
};

//! Set equilibrium model traits
template<class TypeTag>
struct EquilibriumModelTraits<TypeTag, TTag::TwoPTwoCNonEquil>
{
private:
     using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
     using EquilibriumTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = TwoPTwoCUnconstrainedModelTraits<EquilibriumTraits>;
};

//! In case we do not assume full thermal non-equilibrium (e.g. only an energy balance for the solid phase and a fluid mixture)
//! one needs a law for calculating the thermal conductivity of the fluid mixture
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPTwoCNonEquil>
{ using type = ThermalConductivitySimpleFluidLumping<GetPropType<TypeTag, Properties::Scalar>>; };

//! Use the nonequilibrium volume variables together with the 2p2c vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoCNonEquil>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
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

    static constexpr bool useConstraintSolver = getPropValue<TypeTag, Properties::UseConstraintSolver>();
    using EquilibriumVolVars = TwoPTwoCVolumeVariables<NCTraits<BaseTraits, DT, EDM>, useConstraintSolver>;
public:
    using type = NonEquilibriumVolumeVariables<NCTraits<BaseTraits, DT, EDM>, EquilibriumVolVars>;
};

/////////////////////////////////////////////////
// Properties for the non-equilibrium nonisothermal TwoPTwoC model which assumes thermal equilibrium but chemical nonequilibrium
/////////////////////////////////////////////////

namespace TTag {
struct TwoPTwoCNINonEquil { using InheritsFrom = std::tuple<TwoPTwoCNonEquil>; };
} // end namespace TTag

//! Set the non-isothermal model traits with the nonequilibrium model traits as isothermal traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPTwoCNINonEquil>
{
private:
    private:
    using EquiTraits = GetPropType<TypeTag, Properties::EquilibriumModelTraits>;
    static constexpr bool enableTNE = getPropValue<TypeTag, Properties::EnableThermalNonEquilibrium>();
    static constexpr bool enableCNE = getPropValue<TypeTag, Properties::EnableChemicalNonEquilibrium>();
    static constexpr int numEF = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr int numES = getPropValue<TypeTag, Properties::NumEnergyEqSolid>();
    static constexpr auto nf = getPropValue<TypeTag, Properties::NusseltFormulation>();
    static constexpr auto ns = getPropValue<TypeTag, Properties::SherwoodFormulation>();

    using IsothermalTraits = NonEquilibriumModelTraits<EquiTraits, enableCNE, enableTNE, numEF, numES, nf, ns>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the equilibrium IO fields which are in that case the nonisothermal io fields
template<class TypeTag>
struct EquilibriumIOFields<TypeTag, TTag::TwoPTwoCNINonEquil>
{
private:
    using NonisothermalIOFields = EnergyIOFields<TwoPNCIOFields>;
public:
    using type = NonisothermalIOFields;
};

//! Use the nonequilibrium volume variables together with the 2p2c vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPTwoCNINonEquil>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
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
    static constexpr bool useConstraintSolver = getPropValue<TypeTag, Properties::UseConstraintSolver>();
    using EquilibriumVolVars = TwoPTwoCVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>, useConstraintSolver>;
public:
    using type = NonEquilibriumVolumeVariables<NCNITraits<BaseTraits, DT, EDM, ETCM>, EquilibriumVolVars>;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPTwoCNINonEquil>
{ using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Properties
} // end namespace Dumux

#endif
