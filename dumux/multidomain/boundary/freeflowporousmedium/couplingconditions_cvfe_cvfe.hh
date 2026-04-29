// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowPorousMediumCouplingConditions
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_CVFE_CVFE_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_CVFE_CVFE_HH

#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include <dumux/flux/darcyslaw_fwd.hh>
#include <dumux/flux/fickslaw_fwd.hh>
#include <dumux/flux/forchheimerslaw_fwd.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/traits.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/indexhelper.hh>

// ToDo: Needed for coupling options. Maybe export options to seperate file
#include "couplingconditions_staggered_cctpfa.hh"

namespace Dumux {

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class FFPMCouplingConditionsCvfeImpl;

/*!
* \ingroup FreeFlowPorousMediumCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using FFPMCouplingConditionsCvfe
    = FFPMCouplingConditionsCvfeImpl<
        MDTraits, CouplingManager,
        GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
        (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)
    >;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class FFPMCouplingConditionsCvfeImplBase
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FluidSystem = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id> using ModelTraits = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

public:
    static constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;
private:

    using AdvectionType = GetPropType<SubDomainTypeTag<porousMediumIndex>, Properties::AdvectionType>;
    using DarcysLaw = Dumux::DarcysLaw<SubDomainTypeTag<porousMediumIndex>>;
    using ForchheimersLaw = Dumux::ForchheimersLaw<SubDomainTypeTag<porousMediumIndex>>;

    static constexpr bool adapterUsed = ModelTraits<porousMediumIndex>::numFluidPhases() > 1;
    using IndexHelper = FreeFlowPorousMediumCoupling::IndexHelper<freeFlowMassIndex, porousMediumIndex, FluidSystem<freeFlowMassIndex>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<freeFlowMassIndex>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<porousMediumIndex>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(FreeFlowPorousMediumCoupling::IsSameFluidSystem<FluidSystem<freeFlowMassIndex>,
                  FluidSystem<porousMediumIndex>>::value,
                  "All submodels must use the same fluid system");

    using DiffusionCoefficientAveragingType = typename FreeFlowPorousMediumCouplingOptions::DiffusionCoefficientAveragingType;

public:
    /*!
     * \brief Returns the corresponding phase index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i> id, int coupledPhaseIdx = 0)
    { return IndexHelper::couplingPhaseIdx(id, coupledPhaseIdx); }

    /*!
     * \brief Returns the corresponding component index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i> id, int coupledCompIdx)
    { return IndexHelper::couplingCompIdx(id, coupledCompIdx); }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    template<class IpData, class Context>
    static auto darcyPermeability(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                  const IpData& ipData,
                                  const Context& context)
    {
        const auto& pmFvGeometry = context.fvGeometry;
        const auto& pos = ipData.global();

        // TODO: Maybe we find a better solution, which would require to change permeability interfaces
        // to work with ipData.
        Scalar minDistance = 0.0;
        std::size_t closestScvfIdx = 0;
        bool foundClosestScvf = false;

        for (const auto& pmScvf : scvfs(pmFvGeometry))
        {
            if(!pmScvf.boundary())
                continue;

            const auto distance = (pmScvf.ipGlobal() - pos).two_norm();
            if (!foundClosestScvf || distance < minDistance)
            {
                minDistance = distance;
                closestScvfIdx = pmScvf.index();
                foundClosestScvf = true;
            }
        }

        if (!foundClosestScvf)
            DUNE_THROW(Dune::InvalidStateException, "Expected at least one SCVF in porous-medium element geometry");

        return context.vars(pmFvGeometry.scv(pmFvGeometry.scvf(closestScvfIdx).insideScvIdx())).permeability();
    }

    /*!
     * \brief Returns the porous medium density at a given interpolation point
     */
    template<class FVElementGeometry, class ContextOrElementVars, class IpData>
    static Scalar density(const FVElementGeometry& fvGeometry,
                          const ContextOrElementVars& vars,
                          const IpData& ipData)
    {
        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);

        const auto& localBasis = fvGeometry.feLocalBasis();

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        localBasis.evaluateFunction(ipData.local(), shapeValues);

        // interpolate density
        // TODO: Replace by localDof loop, for now we assume that the mass schemes are pure fv schemes
        Scalar density(0.0);
        for (const auto& scv : scvs(fvGeometry))
            density += vars[scv].density(pmPhaseIdx)*shapeValues[scv.localDofIndex()][0];

        return density;
    }

    /*!
     * \brief Returns the porous medium pressure at a given interpolation point
     */
    template<class FVElementGeometry, class ContextOrElementVars, class IpData>
    static Scalar pressure(const FVElementGeometry& fvGeometry,
                           const ContextOrElementVars& vars,
                           const IpData& ipData)
    {
        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);

        const auto& localBasis = fvGeometry.feLocalBasis();

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        localBasis.evaluateFunction(ipData.local(), shapeValues);

        // interpolate pressure
        // TODO: Replace by localDof loop, for now we assume that the mass schemes are pure fv schemes
        Scalar pressure(0.0);
        for (const auto& scv : scvs(fvGeometry))
            pressure += vars[scv].pressure(pmPhaseIdx)*shapeValues[scv.localDofIndex()][0];

        return pressure;
    }

    /*!
     * \brief Evaluate an advective flux across the interface and consider upwinding.
     */
    static Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream)
    {
        const Scalar upwindWeight = 1.0; //TODO use Flux.UpwindWeight or something like Coupling.UpwindWeight?

        if(insideIsUpstream)
            return (upwindWeight * insideQuantity + (1.0 - upwindWeight) * outsideQuantity) * volumeFlow;
        else
            return (upwindWeight * outsideQuantity + (1.0 - upwindWeight) * insideQuantity) * volumeFlow;
    }

protected:

    /*!
     * \brief Returns the distance between an scvf and the corresponding scv center.
     */
    template<class Scv, class Scvf>
    Scalar getDistance_(const Scv& scv, const Scvf& scvf) const
    {
        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
    }
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class FFPMCouplingConditionsCvfeImpl<MDTraits, CouplingManager, enableEnergyBalance, false>
: public FFPMCouplingConditionsCvfeImplBase<MDTraits, CouplingManager>
{
    using ParentType = FFPMCouplingConditionsCvfeImplBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

    static_assert(GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");

    using DiffusionCoefficientAveragingType = typename FreeFlowPorousMediumCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<ParentType::porousMediumIndex>& element,
                                   const FVElementGeometry<ParentType::porousMediumIndex>& fvGeometry,
                                   const ElementVolumeVariables<ParentType::porousMediumIndex>& darcyElemVolVars,
                                   const SubControlVolumeFace<ParentType::porousMediumIndex>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        DUNE_THROW(Dune::NotImplemented, "Energy coupling condition");
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<ParentType::freeFlowMassIndex>& element,
                                   const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                   const ElementVolumeVariables<ParentType::freeFlowMassIndex>& stokesElemVolVars,
                                   const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        DUNE_THROW(Dune::NotImplemented, "Energy coupling condition");
    }

};

} // end namespace Dumux

#endif
