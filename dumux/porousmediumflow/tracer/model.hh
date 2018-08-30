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
 * \ingroup TracerModel
 * \brief Adaption of the fully implicit scheme to the tracer transport model.
 *
 * This model implements a transport of a tracer, where density and viscosity of the
 * fluid phase in which the tracer gets transported are not affected by the tracer.
 * The velocity field is a given spatial parameter.
 * The model is mainly useful for fast computations on given or precomputed
 * velocity fields and thus generally makes sense only in combination with an incompressible
 * one-phase flow velocity field or analytically given / artificial fields. However, reactions
 * between multiple tracer components can be implemented.
 *
 * The transport of the components \f$\kappa \in \{ a, b, c, ... \}\f$ is described by the following equation:
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
 * Note that the tracer model is always considered non-isothermal.
 * The velocity output is fully compatible with the tracer model if you want to write the velocity field to vtk.
*/

#ifndef DUMUX_TRACER_MODEL_HH
#define DUMUX_TRACER_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/discretization/stationaryvelocityfield.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/porousmediumflow/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup TracerModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether mole or mass balances are used
 */
template<int nComp, bool useMol>
struct TracerModelTraits
{
    using Indices = TracerIndices;

    static constexpr int numEq() { return nComp; }
    static constexpr int numPhases() { return 1; }
    static constexpr int numComponents() { return nComp; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool useMoles() { return useMol; }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = useMoles() ? "x" : "X";
        return xString + "^" + FluidSystem::componentName(pvIdx);
    }
};

/*!
 * \ingroup TracerModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam MT The model traits
 */
template<class PV, class FSY, class SSY, class SST, class MT>
struct TracerVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using SolidSystem = SSY;
    using SolidState = SST;
    using ModelTraits = MT;
};

// \{
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit tracer model.
NEW_TYPE_TAG(Tracer, INHERITS_FROM(PorousMediumFlow));

///////////////////////////////////////////////////////////////////////////
// properties for the tracer model
///////////////////////////////////////////////////////////////////////////

//! Define that mole fractions are used in the balance equations
SET_BOOL_PROP(Tracer, UseMoles, true);

//! set the model traits
SET_PROP(Tracer, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = TracerModelTraits<FluidSystem::numComponents, GET_PROP_VALUE(TypeTag, UseMoles)>;
};

//! Use the tracer local residual function for the tracer model
SET_TYPE_PROP(Tracer, LocalResidual, TracerLocalResidual<TypeTag>);

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(Tracer, VtkOutputFields, TracerVtkOutputFields);

//! Set the volume variables property
SET_PROP(Tracer, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = TracerVolumeVariablesTraits<PV, FSY, SSY, SST, MT>;
public:
    using type = TracerVolumeVariables<Traits>;
};

//! We use darcy's law as the default for the advective fluxes
SET_TYPE_PROP(Tracer, AdvectionType, StationaryVelocityField<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Use simple model with constant tortuosity as pm diffusivity model
SET_TYPE_PROP(Tracer, EffectiveDiffusivityModel, DiffusivityConstantTortuosity<typename GET_PROP_TYPE(TypeTag, Scalar)>);
} // end namespace Properties
// \}
} // end namespace Dumux

#endif
