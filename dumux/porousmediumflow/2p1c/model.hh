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
 * \ingroup TwoPOneCModel
 * \brief A two-phase one-component flow model using the fully implicit scheme.
 *
 * \note The 2p1c model requires the use of the non-isothermal extension found in dumux/porousmediumflow/nonisothermal.
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
 * If only one phase is present, \f$p_n\f$ and \f$T\f$ are the primary variables.
 * In the presence of two phases, \f$p_n\f$ and \f$S_w\f$ become the new primary variables.
 */

#ifndef DUMUX_2P1C_MODEL_HH
#define DUMUX_2P1C_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "darcyslaw.hh"
#include "iofields.hh"
#include "localresidual.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "primaryvariableswitch.hh"

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief Specifies a number properties of models
 *        considering two phases with water as a single component.
 */
template<TwoPFormulation f>
struct TwoPOneCNIModelTraits
{
    using Indices = TwoPOneCIndices;

    //! We solve for one more equation, i.e. the energy balance
    static constexpr int numEq() { return 2; }
    //! only one energy equation is needed when assuming thermal equilibrium
    static constexpr int numEnergyEq() { return 1; }
    static constexpr int numPhases() { return 2; }
    static constexpr int numComponents() { return 1; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return true; }

    static constexpr TwoPFormulation priVarFormulation() { return f; }

    template <class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        if (priVarFormulation() == TwoPFormulation::p0s1)
            return (pvIdx == 0) ? "p_w" :
                                  (state == Indices::twoPhases) ? "S_n" : "T";
        else
            return (pvIdx == 0) ? "p_n" :
                                  (state == Indices::twoPhases) ? "S_w" : "T";
    }
};

/*!
 * \ingroup TwoPOneCModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct TwoPOneCNIVolumeVariablesTraits
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
//! The type tag for the non-isothermal two-phase one-component model.
NEW_TYPE_TAG(TwoPOneCNI, INHERITS_FROM(PorousMediumFlow));

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

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
};

//! Set the default formulation to pw-sn
SET_PROP(TwoPOneCNI, Formulation)
{ static constexpr TwoPFormulation value = TwoPFormulation::p1s0; };

//! Do not block spurious flows by default.
SET_BOOL_PROP(TwoPOneCNI, UseBlockingOfSpuriousFlow, false);

//! The specific local residual (i.e. balance equations).
SET_TYPE_PROP(TwoPOneCNI, LocalResidual, TwoPOneCLocalResidual<TypeTag>);

//! Use a modified version of Darcy's law which allows for blocking of spurious flows.
SET_TYPE_PROP(TwoPOneCNI, AdvectionType, TwoPOneCDarcysLaw<TypeTag>);

//! Set the volume variables property
SET_PROP(TwoPOneCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    static_assert(FSY::numComponents == 1, "Only fluid systems with 1 component are supported by the 2p1cni model!");
    static_assert(FSY::numPhases == 2, "Only fluid systems with 2 phases are supported by the 2p1cni model!");

    using Traits = TwoPOneCNIVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = TwoPOneCVolumeVariables<Traits>;
};

//! The primary variable switch for the 2p1cni model.
SET_TYPE_PROP(TwoPOneCNI, PrimaryVariableSwitch, TwoPOneCPrimaryVariableSwitch);

//! The primary variables vector for the 2p1cni model.
SET_PROP(TwoPOneCNI, PrimaryVariables)
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                     GET_PROP_TYPE(TypeTag, ModelTraits)::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity.
SET_TYPE_PROP(TwoPOneCNI, ThermalConductivityModel, ThermalConductivitySomerton<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//! Set the non-isothermal model traits
SET_TYPE_PROP(TwoPOneCNI, ModelTraits, TwoPOneCNIModelTraits<GET_PROP_VALUE(TypeTag, Formulation)>);

//! The non-isothermal vtk output fields.
SET_PROP(TwoPOneCNI, IOFields)
{
    using TwoPOneCIOF = TwoPOneCIOFields<GET_PROP_VALUE(TypeTag, Formulation), TwoPOneCIndices>;
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using type = EnergyIOFields<TwoPOneCIOF, ModelTraits>;
};

} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_2P1C_MODEL_HH
