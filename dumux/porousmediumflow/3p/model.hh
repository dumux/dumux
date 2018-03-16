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
#include <dumux/porousmediumflow/nonisothermal/vtkoutputfields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Specifies a number properties of three-phase models.
 */
struct ThreePModelTraits
{
    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 3; }
    static constexpr int numComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
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

//! The VolumeVariables property
SET_TYPE_PROP(ThreeP, VolumeVariables, ThreePVolumeVariables<TypeTag>);

//! The indices required by the isothermal 3p model
SET_PROP(ThreeP, Indices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ThreePIndices<FluidSystem,/*PVOffset=*/0>;
};

//! The spatial parameters to be employed.
//! Use FVSpatialParams by default.
SET_TYPE_PROP(ThreeP, SpatialParams, FVSpatialParams<TypeTag>);

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
SET_PROP(ThreeP, VtkOutputFields)
{
private:
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = ThreePVtkOutputFields<Indices>;
};

/////////////////////////////////////////////////
// Properties for the non-isothermal 3p model
/////////////////////////////////////////////////

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySomerton<Scalar, Indices>;
};

//! Set non-isothermal output fields
SET_PROP(ThreePNI, VtkOutputFields)
{
private:
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using IsothermalFields = ThreePVtkOutputFields<Indices>;
public:
    using type = EnergyVtkOutputFields<IsothermalFields>;
};

//! Set non-isothermal Indices
SET_PROP(ThreePNI, Indices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using IsothermalIndices = ThreePIndices<FluidSystem,/*PVOffset=*/0>;
    static constexpr int numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
public:
    using type = EnergyIndices<IsothermalIndices, numEq, 0>;
};

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
