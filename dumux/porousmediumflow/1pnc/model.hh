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
 * \ingroup OnePNCModel
 * \brief  Adaption of the fully implicit model to the one-phase n-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists
 * of n components, using a standard Darcy approach as the equation for the
 * conservation of momentum:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \phi\frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a, ... \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p -
 \varrho {\textbf g} \right)
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fraction of dissolved components \f$x^\kappa\f$.
 */

#ifndef DUMUX_1PNC_MODEL_HH
#define DUMUX_1PNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/nonequilibrium/model.hh>
#include <dumux/porousmediumflow/nonequilibrium/volumevariables.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Specifies a number properties of models that
 *        consider a single-phase with multiple components.
 *
 * \tparam nComp the number of components to be considered.
 */
template<int nComp, bool useM, int repCompEqIdx = nComp>
struct OnePNCModelTraits
{
    using Indices = OnePNCIndices;

    static constexpr int numEq() { return nComp; }
    static constexpr int numPhases() { return 1; }
    static constexpr int numComponents() { return nComp; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    static constexpr bool useMoles() { return useM; }
    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = useMoles() ? "x" : "X";
        if (pvIdx == 0)
            return "p";
        else
            return xString + "^" + FluidSystem::componentName(pvIdx)
                   + "_" + FluidSystem::phaseName(0);
    }
};

/*!
 * \ingroup OnePNCModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct OnePNCVolumeVariablesTraits
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

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the implicit the isothermal & non-isothermal one phase n component problems
// Create new type tags
namespace TTag {
struct OnePNC { using InheritsFrom = std::tuple<PorousMediumFlow>; };
struct OnePNCNI { using InheritsFrom = std::tuple<OnePNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! Set as default that no component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::OnePNC> { static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::numComponents; };

//! The base model traits. Per default, we use the number of components of the fluid system.
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::OnePNC>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = OnePNCModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNC> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::OnePNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::OnePNC>
{ using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! Use mole fractions in the balance equations by default
template<class TypeTag>
struct UseMoles<TypeTag, TTag::OnePNC> { static constexpr bool value = true; };

//! The local residual function
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePNC> { using type = CompositionalLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = OnePNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = OnePNCVolumeVariables<Traits>;
};

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::OnePNC> { using type = OnePNCIOFields; };

///////////////////////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! the non-isothermal vtk output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::OnePNCNI> { using type = EnergyIOFields<OnePNCIOFields>; };

//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::OnePNCNI>
{ using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>; };

//! model traits of the non-isothermal model.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNCNI>
{
private:
    using IsothermalTraits = GetPropType<TypeTag, Properties::BaseModelTraits>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

} // end namespace Properties

template<class OnePNCModelTraits>
struct OnePNCNonequilibriumModelTraits : public OnePNCModelTraits
{
    static constexpr int numConstraintEq() { return 0; }
};

namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
namespace TTag {
struct OnePNCNonEquil { using InheritsFrom = std::tuple<NonEquilibrium, OnePNC>; };
} // end namespace TTag

/////////////////////////////////////////////////
// Properties for the non-equilibrium OnePNC model
/////////////////////////////////////////////////

template<class TypeTag>
struct EquilibriumLocalResidual<TypeTag, TTag::OnePNCNonEquil> { using type = CompositionalLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct EquilibriumIOFields<TypeTag, TTag::OnePNCNonEquil> { using type = OnePNCIOFields; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OnePNCNonEquil>
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
// by default chemical non equilibrium is enabled in the nonequil model, switch that off here
template<class TypeTag>
struct EnableChemicalNonEquilibrium<TypeTag, TTag::OnePNCNonEquil> { static constexpr bool value = false; };

//! set equilibrium model traits
template<class TypeTag>
struct EquilibriumModelTraits<TypeTag, TTag::OnePNCNonEquil>
{
private:
     using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
     using EquilibriumTraits = OnePNCModelTraits<FluidSystem::numComponents,  getPropValue<TypeTag, Properties::UseMoles>(), getPropValue<TypeTag, Properties::ReplaceCompEqIdx>()>;
public:
    using type = OnePNCNonequilibriumModelTraits<EquilibriumTraits>;
};

//! in case we do not assume full non-equilibrium one needs a thermal conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::OnePNCNonEquil>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivityAverage<Scalar>;
};

//! use the mineralization volume variables together with the 2pnc vol vars
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNCNonEquil>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    using Traits = OnePNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
    using EquilibriumVolVars = OnePNCVolumeVariables<Traits>;
public:
    using type = NonEquilibriumVolumeVariables<Traits, EquilibriumVolVars>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
