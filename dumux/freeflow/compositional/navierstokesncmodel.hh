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
 * \ingroup FreeflowNCModel
 *
 * \brief A single-phase, multi-component free-flow model
 *
 * For an equations not specific to multiple components see dumux/freeflow/navierstokes/model.hh
 *
 * The multi-component system is closed by a <B> component mass/mole balance equation </B> for each component \f$\kappa\f$:
 * \f[
 *    \frac{\partial \left(\varrho X^\kappa\right)}{\partial t}
 *    + \nabla \cdot \left( \varrho {\boldsymbol{v}} X^\kappa
 *    - (D^\kappa + D_\text{t}) \varrho \textbf{grad}\, X^\kappa \right)
 *    - q^\kappa = 0
 * \f]
 *
 * Alternatively, one component balance equation can be replace by a <B> total mass/mole balance equation </B>:
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \nabla \cdot \left(
 *        \varrho {\boldsymbol{v}}
 *        - \sum_\kappa (D^\kappa + D_\text{t}) \varrho \textbf{grad}\, X^\kappa
 *      \right)
 *    - q = 0
 * \f]
 *
 * The eddy diffusivity \f$ D_\text{t} \f$ is related to the eddy viscosity \f$ \nu_\text{t} \f$
 * by the turbulent Schmidt number, for Navier-Stokes models \f$ D_\text{t} = 0 \f$.
 * \f[ D_\text{t} = \frac{\nu_\text{t}}{\mathrm{Sc}_\text{t}} \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_FREEFLOW_NC_MODEL_HH
#define DUMUX_FREEFLOW_NC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/freeflow/nonisothermal/iofields.hh>
#include <dumux/flux/fickslaw.hh>
#include <dumux/flux/fourierslaw.hh>

#include "volumevariables.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "iofields.hh"

#include <dumux/assembly/staggeredlocalresidual.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the multi-component free-flow model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam useM Use molar or mass balances
 * \tparam repCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 */
template<int dimension, int nComp, bool useM, int repCompEqIdx = nComp>
struct NavierStokesNCModelTraits : NavierStokesModelTraits<dimension>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp; }

    //! The number of components
    static constexpr int numFluidComponents() { return nComp; }

    //! Use moles or not
    static constexpr bool useMoles() { return useM; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return true; }

    //! Index of of a component balance eq. to be replaced by a total mass/mole balance
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return false; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::none; }

    //! the indices
    using Indices = NavierStokesIndices<dimension>;
};

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component free-flow model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, multi-component isothermal free-flow model
struct NavierStokesNC { using InheritsFrom = std::tuple<FreeFlow>; };

//! The type tag for the single-phase, multi-component non-isothermal free-flow model
struct NavierStokesNCNI { using InheritsFrom = std::tuple<NavierStokesNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

//!< states some specifics of the free-flow model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesNC>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();

public:
    using type = NavierStokesNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::NavierStokesNC> { static constexpr bool value = false; }; //!< Defines whether molar (true) or mass (false) density is used
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::NavierStokesNC> { static constexpr int value = 0; }; //<! Set the ReplaceCompEqIdx to 0 by default
template<class TypeTag>
struct NormalizePressure<TypeTag, TTag::NavierStokesNC> { static constexpr bool value = true; }; //!< Normalize the pressure term in the momentum balance by default

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesNC> { using type = FreeflowNCResidual<TypeTag>; };

//! Use Fick's law for molecular diffusion per default
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::NavierStokesNC> { using type = FicksLaw<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using BaseTraits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    template<class BaseTraits, class DT>
    struct NCTraits : public BaseTraits { using DiffusionType = DT; };
public:
    using type = FreeflowNCVolumeVariables<NCTraits<BaseTraits, DT>>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesNC> { using type = FreeflowNCFluxVariables<TypeTag>; };

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesNC> { using type = FreeflowNCIOFields<NavierStokesIOFields>; };

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::NavierStokesNC>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component free-flow model
//////////////////////////////////////////////////////////////////////////

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesNCNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
    using IsothermalModelTraits = NavierStokesNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! The non-isothermal I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesNCNI>
{
private:
    using IsothermalIOFields = FreeflowNCIOFields<NavierStokesIOFields>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalIOFields>;
};

//! Use Fourier's Law as default heat conduction type
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesNCNI> { using type = FouriersLaw<TypeTag>; };

} // end namespace Properties
} // end namespace Dumux

#endif
