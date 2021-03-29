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
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_NAVIERSTOKES_1P_MODEL_HH
#define DUMUX_NAVIERSTOKES_1P_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "./../../iofields.hh"

#include <dumux/flux/cctpfa/fourierslaw.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/energy/model.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the Navier-Stokes model
 *
 * \tparam dimension The dimension of the problem
 */
struct NavierStokesMassOnePModelTraits
{
    //! There are as many momentum balance equations as dimensions
    //! and one mass balance equation.
    //! plus one equation per phasefield and per transported species
    static constexpr int numEq() { return 1+3+3; }

    //! The number of phases is 1
    static constexpr int numFluidPhases() { return 1; }

    //! The number of components is 1
    static constexpr int numFluidComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return false; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return false; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::none; }

    //! the indices
    using Indices = NavierStokesMassOnePIndices;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the volume variables of the Navier-Stokes model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class FST,
         class MT>
struct NavierStokesMassOnePVolumeVariablesTraits
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
struct NavierStokesMassOneP{ using InheritsFrom = std::tuple<ModelProperties>; };
struct NavierStokesMassOnePNI{ using InheritsFrom = std::tuple<NavierStokesMassOneP>; };

}


//!< states some specifics of the Navier-Stokes model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneP>
{ using type = NavierStokesMassOnePModelTraits; };

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::NavierStokesMassOneP>{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOneP>
{ using type = NavierStokesMassOnePLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMassOnePVolumeVariables<Traits>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::NavierStokesMassOneP> { using type = NavierStokesMassOnePFluxVariables<TypeTag>; };

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOneP> { using type = FluxVariablesCaching::EmptyCache<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOneP> { using type = FluxVariablesCaching::EmptyCacheFiller; };

// ! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneP> { using type = NavierStokesIOFields; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::NavierStokesMassOneP>
{
private:
    struct EmptyCouplingManager {};
public:
    using type = EmptyCouplingManager;
};

///////////////////////////////////////////////////////////////////////////
// Properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! Add temperature to the output
template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOnePNI> { using type = NavierStokesEnergyIOFields<NavierStokesIOFields>; };

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOnePNI> { using type = NavierStokesEnergyModelTraits<NavierStokesMassOnePModelTraits>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOnePNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using BaseTraits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
    using ETCM = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };
public:
    using type = NavierStokesMassOnePVolumeVariables<NITraits>;
};

//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOnePNI>
{
    struct type
    {
        template<class VolVars>
        static auto effectiveThermalConductivity(const VolVars& volVars)
        {
            return volVars.fluidThermalConductivity();
        }
    };
    // using type = ThermalConductivityAverage<GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NavierStokesMassOnePNI>
{ using type = FouriersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>; };

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::NavierStokesMassOnePNI>
{
    struct Cache : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache
    {

    };

    using type = Cache;
};

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::NavierStokesMassOnePNI>
{
    class Filler : public GetPropType<TypeTag, Properties::HeatConductionType>::Cache::Filler
    {
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
        using GridView = typename GridGeometry::GridView;
        using Element = typename GridView::template Codim<0>::Entity;

        using FVElementGeometry = typename GridGeometry::LocalView;
        using SubControlVolume = typename GridGeometry::SubControlVolume;
        using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
        using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    public:
        static constexpr bool isSolDependent = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

        Filler(const Problem& problem)
        : problemPtr_(&problem) {}

        template<class FluxVariablesCacheContainer, class FluxVariablesCache>
        void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                  FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf,
                  bool forceUpdateAll = false)
        {
            using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
            using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

            // forward to the filler of the diffusive quantities
            HeatConductionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
        }

    private:
        const Problem& problem() const
        { return *problemPtr_; }

        const Problem* problemPtr_;
    };

    using type = Filler;
};

template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::NavierStokesMassOnePNI> { static constexpr bool value = true; };

} // end namespace Properties
// }

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_MODEL_HH
