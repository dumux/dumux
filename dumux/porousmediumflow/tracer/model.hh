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
 + \varrho D^\kappa_\text{pm} \textbf{grad} X^\kappa \right\rbrace = q.
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
#include <dumux/flux/stationaryvelocityfield.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/porousmediumflow/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
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
    static constexpr int numFluidPhases() { return 1; }
    static constexpr int numFluidComponents() { return nComp; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool useMoles() { return useMol; }
};

/*!
 * \ingroup TracerModel
 * \brief Traits class for the volume variables of the single-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam MT The model traits
 */
template<class PV, class FSY, class SSY, class SST, class MT, class DT, class EDM>
struct TracerVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using SolidSystem = SSY;
    using SolidState = SST;
    using ModelTraits = MT;
    using DiffusionType = DT;
    using EffectiveDiffusivityModel = EDM;
};

// \{
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit tracer model.
// Create new type tags
namespace TTag {
struct Tracer { using InheritsFrom = std::tuple<PorousMediumFlow>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the tracer model
///////////////////////////////////////////////////////////////////////////

//! Define that mole fractions are used in the balance equations
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Tracer> { static constexpr bool value = true; };

//! set the model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Tracer>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = TracerModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>()>;
};

//! Use the tracer local residual function for the tracer model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Tracer> { using type = TracerLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::Tracer> { using type = TracerIOFields; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Tracer>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

    using Traits = TracerVolumeVariablesTraits<PV, FSY, SSY, SST, MT, DT, EDM>;
public:
    using type = TracerVolumeVariables<Traits>;
};

//! We use darcy's law as the default for the advective fluxes
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::Tracer> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };

//! Use simple model with constant tortuosity as pm diffusivity model
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::Tracer> { using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };
} // end namespace Properties
// \}
} // end namespace Dumux

#endif
