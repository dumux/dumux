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
 * \ingroup AdvectionDiffusionModel
 * \brief The advection-diffusion transport Model
 *
 * The transport balance is described by the following equation:
 \f[
    \frac{\partial u}{\partial t} +
    \nabla \cdot (u \textbf{v}) -
    \nabla \cdot (D \nabla u) = q
 \f]
 * where \f$u\f$ is a transported scalar,
 * \f$\textbf{v}\f$ is a constant velocity field,
 * and \f$D\f$ is a constant isometric diffusive coefficient.
*/

#ifndef DUMUX_ADVECTION_DIFFUSION_MODEL_HH
#define DUMUX_ADVECTION_DIFFUSION_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/flux/stationaryvelocityfield.hh>
#include <dumux/flux/fickslaw.hh>

#include "iofields.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup AdvectionDiffusionModel
 * \brief The Indices
 */
struct AdvectionDiffusionIndices
{
    static constexpr int componentIdx = 0;
    static constexpr int transportEqIdx = componentIdx;
};

/*!
 * \ingroup AdvectionDiffusionModel
 * \brief The Model Traits
 */
template<int nComp, bool useMol>
struct AdvectionDiffusionModelTraits
{
    using Indices = AdvectionDiffusionIndices;

    static constexpr int numFluidComponents() { return nComp; }
    static constexpr int numEq() { return nComp; }
    static constexpr int numFluidPhases() { return 1; }
    static constexpr bool useMoles() { return useMol; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup AdvectionDiffusionModel
 * \brief The volume variable traits
 */
template<class PV, class FSY, class FST, class MT, class DT, class AT>
struct AdvectionDiffusionVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using ModelTraits = MT;
    using DiffusionType = DT;
    using AdvectionType = AT;
};

// \{
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit tracer model.
// Create new type tags
namespace TTag {
struct AdvectionDiffusion { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// Properties
///////////////////////////////////////////////////////////////////////////

//! Use mole fractions in the balance equations by default
template<class TypeTag>
struct UseMoles<TypeTag, TTag::AdvectionDiffusion> { static constexpr bool value = true; };

//!< states some specifics of the free-flow model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::AdvectionDiffusion>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

public:
    using type = AdvectionDiffusionModelTraits<numComponents, useMoles>;
};

//! Use Fick's law for molecular diffusion per default
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::AdvectionDiffusion> { using type = FicksLaw<TypeTag>; };

//! We use darcy's law as the default for the advective fluxes
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::AdvectionDiffusion> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };

//! set the local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::AdvectionDiffusion>
{ using type = AdvectionDiffusionLocalResidual<TypeTag>;};

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::AdvectionDiffusion> { using type = AdvectionDiffusionFluxVariables<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::AdvectionDiffusion> { using type = AdvectionDiffusionIOFields; };

template<class TypeTag>
struct FluidState<TypeTag, TTag::AdvectionDiffusion>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::AdvectionDiffusion>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using AT = GetPropType<TypeTag, Properties::AdvectionType>;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using Traits = AdvectionDiffusionVolumeVariablesTraits<PV, FSY, FST, MT, DT, AT>;

public:
    using type = AdvectionDiffusionVolumeVariables<Traits>;
};

} // end namespace Properties
// \}
} // end namespace Dumux

#endif
