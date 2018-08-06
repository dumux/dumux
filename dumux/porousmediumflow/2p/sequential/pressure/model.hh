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
#ifndef DUMUX_PRESSURE_MODEL_HH
#define DUMUX_PRESSURE_MODEL_HH

/*!
 * \file
 * \ingroup PressureModel
 * \brief Adaption of the fully implicit scheme to the Pressure Pressure model.
 *
 * This model implements a Pressure of a Pressure, where density and viscosity of the
 * fluid phase in which the Pressure gets Pressureed are not affected by the Pressure.
 * The velocity field is a given spatial parameter.
 * The model is mainly useful for fast computations on given or precomputed
 * velocity fields and thus generally makes sense only in combination with an incompressible
 * one-phase flow velocity field or analytically given / artificial fields. However, reactions
 * between multiple Pressure components can be implemented.
 *
 * The Pressure of the components \f$\kappa \in \{ a, b, c, ... \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa {\textbf v_f}
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme (TPFA or MPFA) as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables the mole or mass fraction of dissolved components \f$x\f$.
 * Note that the Pressure model is always considered non-isothermal.
 * The velocity output is fully compatible with the Pressure model if you want to write the velocity field to vtk.
*/

#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/discretization/stationaryvelocityfield.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/porousmediumflow/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup PressureModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether mole or mass balances are used
 */
struct PressureModelTraits
{
    using Indices = PressureIndices;

    static constexpr int numEq() { return 1; }
    static constexpr int numPhases() { return 2; }
    static constexpr int numComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup PressureModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT>
struct PressureVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

// \{
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit Pressure model.
NEW_TYPE_TAG(Pressure, INHERITS_FROM(PorousMediumFlow));

///////////////////////////////////////////////////////////////////////////
// properties for the Pressure model
///////////////////////////////////////////////////////////////////////////

//! set the model traits
SET_PROP(Pressure, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = PressureModelTraits;
};

//! Use the Pressure local residual function for the Pressure model
SET_TYPE_PROP(Pressure, LocalResidual, PressureLocalResidual<TypeTag>);

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(Pressure, VtkOutputFields, PressureVtkOutputFields);

//! Set the volume variables property
SET_PROP(Pressure, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    using Traits = PressureVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = PressureVolumeVariables<Traits>;
};

//! The two-phase model uses the immiscible fluid state
SET_PROP(Pressure, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

} // end namespace Properties
// \}
} // end namespace Dumux

#endif
