// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 *
 * \brief A two-dimensional shallow water equations model
 *
 * The two-dimensional shallow water equations (SWEs) can be written as:
 *
 * \f[
 * \frac{\partial \mathbf{U}}{\partial t} +
 * \frac{\partial \mathbf{F}}{\partial x} +
 * \frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
 * \f]
 *
 * The first equation is the water balance equation (volume balance) and the following two equations balance the
 * momentum in x-direction and y-direction. \f$ \mathbf{U} \f$, \f$ \mathbf{F} \f$ and \f$ \mathbf{G} \f$ are defined as:
 *
 * \f[
 * \mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
 * \mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 - \nu\frac{\partial uh}{\partial x} \\ huv - \nu\frac{\partial vh}{\partial x} \end{bmatrix},
 * \mathbf{G} = \begin{bmatrix} hv \\ huv - \nu\frac{\partial uh}{\partial y} \\ hv^2  + \frac{1}{2} gh^2 - \nu\frac{\partial vh}{\partial y} \end{bmatrix}
 * \f]
 *
 * \f$ h \f$ is the water depth (in \f$ m \f$), \f$ u \f$ the velocity in
 * x-direction and \f$ v \f$ the velocity in y-direction (in \f$ ms^{-1} \f$). \f$ g \f$ is the constant of gravity (in \f$ ms^{-2} \f$).
 * \f$ \nu \f$ is the effective turbulent viscosity (in \f$ m^2s^{-1} \f$).
 * By default the shallow water model neglects the viscous terms, but they can be enabled by setting the parameter `ShallowWater.EnableViscousFlux = true`.
 *
 * The source terms for bed slope \f$ \mathbf{S_b} \f$ and bottom friction \f$ \mathbf{S_f} \f$
 * are given as:
 * \f[
 * \mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
 *                \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
 * \mathbf{S_f} = \begin{bmatrix} 0 \\ -\frac{\tau_{x}}{\rho} \\ -\frac{\tau_{y}}{\rho}\end{bmatrix}.
 * \f]
 *
 * \f$ z \f$ is the bed surface (in \f$ m \f$), \f$ \tau_{x} \f$ and \f$ \tau_{y} \f$  the bottom shear stress (in \f$ Nm^{-2} \f$) in x- an y-direction, respectively.
 * \f$ \rho \f$ is the water density (in \f$ kgm^{-3} \f$). The bed slope source term \f$ \mathbf{S_b} \f$ is covered by the hydrostatic reconstruction
 * within the flux computation and must therefore not be treated separately within the computation of the source terms.
 * On the contrary, the bottom friction source term must be implemented in the problem (have a look at <a href="https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/shallowwaterfriction/README.md" target="_blank">the shallow water example</a>).
 *
 * When using a fully implicit Euler time discretization, note the following:
 * Large time step sizes (CFL > 1) can lead to a strong smearing of sharp fronts.
 * This can be seen in the movement of fast traveling waves (e.g. dam break
 * waves). Nevertheless, the fully implicit time discretization shows
 * good results in cases where slowly moving waves are considered. Thus, the model
 * is a good choice for simulating flow in rivers and channels, where the
 * fully-implicit discretization allows large time steps and reduces the
 * overall computation time drastically.
 *
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_MODEL_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/flux/shallowwaterflux.hh>
#include <dumux/flux/shallowwaterviscousflux.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterModels
 * \brief Specifies a number properties of shallow water models.
 */
struct ShallowWaterModelTraits
{
    using Indices = ShallowWaterIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }
};

/*!
 * \ingroup ShallowWaterModels
 * \brief Traits class for the volume variables of the shallow water model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class MT>
struct ShallowWaterVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using ModelTraits = MT;
};


namespace Properties {

//! Type tag for shallow water equation model inherits from model properties
namespace TTag {
struct ShallowWater { using InheritsFrom = std::tuple<ModelProperties>; };
}// end namespace TTag

//////////////////////////////////////////////////////////////////
// Define properties
//////////////////////////////////////////////////////////////////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterResidual<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterFluxVariables<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ShallowWater>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = ShallowWaterVolumeVariablesTraits<PV, FSY, MT>;
public:
    using type = ShallowWaterVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterIOFields; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterFlux< Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>> >; };

template<class TypeTag>
struct ViscousFluxType<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterViscousFlux< Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>> >; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ShallowWater>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

} // end properties
} // end namespace Dumux

#endif
