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

#include <array>

// property forward declarations
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Specifies a number properties of two-phase two-component models.
 *
 * \tparam f The two-phase formulation used
 * \tparam useM Boolean to specify if moles or masses are balanced
 * \tparam replCompEqIdx the equation which is replaced by the total mass balance (none if replCompEqIdx >= numComponents)
 */
template<TwoPFormulation formulation, bool useMol, int replCompEqIdx = 2>
struct TwoPTwoCModelTraits : public TwoPNCModelTraits</*numComps=*/2, useMol, /*setMFracForFirstPhase=*/true, formulation, replCompEqIdx>
{
    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        static const std::string xString = useMol ? "x" : "X";
        static const std::array<std::string, 3> p0s1SwitchedPvNames = {{
            xString + "^" + FluidSystem::componentName(1) + "_" + FluidSystem::phaseName(0),
            xString + "^" + FluidSystem::componentName(0) + "_" + FluidSystem::phaseName(1),
            "S_n"}};
        static const std::array<std::string, 3> p1s0SwitchedPvNames = {{
            xString + "^" + FluidSystem::componentName(1) + "_" + FluidSystem::phaseName(0),
            xString + "^" + FluidSystem::componentName(0) + "_" + FluidSystem::phaseName(1),
            "S_w"}};

        switch (formulation)
        {
        case TwoPFormulation::p0s1:
            return pvIdx == 0 ? "p_w" : p0s1SwitchedPvNames[state-1];
        case TwoPFormulation::p1s0:
            return pvIdx == 0 ? "p_n" : p1s0SwitchedPvNames[state-1];
        }
    }
};

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
    using type = TwoPTwoCModelTraits< getPropValue<TypeTag, Properties::Formulation>(),
                                      getPropValue<TypeTag, Properties::UseMoles>(),
                                      getPropValue<TypeTag, Properties::ReplaceCompEqIdx>() >;
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
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    static_assert(FSY::numComponents == 2, "Only fluid systems with 2 components are supported by the 2p2c model!");
    static_assert(FSY::numPhases == 2, "Only fluid systems with 2 phases are supported by the 2p2c model!");

    static constexpr bool useConstraintSolver = getPropValue<TypeTag, Properties::UseConstraintSolver>();

    using Traits = TwoPNCVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = TwoPTwoCVolumeVariables<Traits, useConstraintSolver>;
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

//! Set non-isothermal output fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPTwoCNI> { using type = EnergyIOFields<TwoPNCIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPTwoCNI> { using type = ThermalConductivitySomerton<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Properties
} // end namespace Dumux

#endif
