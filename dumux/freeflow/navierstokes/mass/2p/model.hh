// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the runtime parameter <code>Problem.EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}})
 * \f]
 * is neglected.
 *
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 */

#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_2P_MODEL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_2P_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariablescachefiller.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the two-phase flow Navier-Stokes mass model
 */
struct NavierStokesMassTwoPModelTraits
{
    //! There are as many momentum balance equations as dimensions
    //! and one mass balance equation.
    static constexpr int numEq() { return 3; }

    //! The number of phases is 1
    static constexpr int numFluidPhases() { return 1; }

    //! The number of components is 1
    static constexpr int numFluidComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The two-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return false; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! the indices
    using Indices = NavierStokesMassTwoPIndices;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the volume variables of the Navier-Stokes model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class FST,
         class MT>
struct NavierStokesMassTwoPVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using ModelTraits = MT;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Navier-Stokes model
struct NavierStokesMassTwoP { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassTwoP>
{ using type = NavierStokesMassTwoPModelTraits; };

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassTwoP>
{ using type = NavierStokesMassTwoPLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassTwoP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");
    using Traits = NavierStokesMassTwoPVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMassTwoPVolumeVariables<Traits>;
};

// The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassTwoP> { using type = NavierStokesMassTwoPIOFields; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::NavierStokesMassTwoP>
{
public:
    using type = struct EmptyCouplingManager {};
};

// Set the default spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::NavierStokesMassTwoP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FreeFlowDefaultSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
