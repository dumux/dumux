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
 * \ingroup ThreePModel
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * This model implements three-phase flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$.
 * The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equations for the conservation
 * of the phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are gas phase pressure \f$p_g\f$,
 * water saturation \f$S_w\f$ and NAPL saturation \f$S_n\f$.
 */
#ifndef DUMUX_3P_MODEL_HH
#define DUMUX_3P_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Specifies a number properties of three-phase models.
 */
struct ThreePModelTraits
{
    using Indices = ThreePIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 3; }
    static constexpr int numComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }

    template <class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        switch (pvIdx)
        {
            case 0: return "p_g";
            case 1: return "S_w";
            default: return "S_n";
        }
    }
};

/*!
 * \ingroup ThreePModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct ThreePVolumeVariablesTraits
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
//! The type tags for the isothermal three-phase model
NEW_TYPE_TAG(ThreeP, INHERITS_FROM(PorousMediumFlow));
//! The type tags for the non-isothermal three-phase model
NEW_TYPE_TAG(ThreePNI, INHERITS_FROM(ThreeP));

//////////////////////////////////////////////////////////////////
// Properties for the isothermal 3p model
//////////////////////////////////////////////////////////////////

//! Set the model traits
SET_PROP(ThreeP, ModelTraits)
{
 private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p model!");
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p model!");
 public:
    using type = ThreePModelTraits;
};

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreeP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);

//! Set the volume variables property
SET_PROP(ThreeP, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    using Traits = ThreePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = ThreePVolumeVariables<Traits>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state.
 *
 *  The fluid state should be chosen appropriately for the model ((non-)isothermal, equilibrium, ...).
 *  This can be done in the problem.
 */
SET_PROP(ThreeP, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(ThreeP, IOFields, ThreePIOFields);

/////////////////////////////////////////////////
// Properties for the non-isothermal 3p model
/////////////////////////////////////////////////

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

//! Set non-isothermal output fields
SET_TYPE_PROP(ThreePNI, IOFields, EnergyIOFields<ThreePIOFields>);

//! Set non-isothermal model traits
SET_PROP(ThreePNI, ModelTraits)
{
private:
   using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
   static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p model!");
   static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p model!");
public:
   using type = PorousMediumFlowNIModelTraits<ThreePModelTraits>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
