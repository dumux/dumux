// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \copydoc Dumux::FreeFlowPoreNetworkCouplingConditions
 */

#ifndef DUMUX_MD_FREEFLOW_PORENETWORK_COUPLINGCONDITIONS_HH
#define DUMUX_MD_FREEFLOW_PORENETWORK_COUPLINGCONDITIONS_HH

#include <numeric>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/fickslaw_fwd.hh>
#include <dumux/freeflow/navierstokes/momentum/velocitygradients.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/traits.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/indexhelper.hh>

namespace Dumux {

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class FreeFlowPoreNetworkCouplingConditionsImplementation;

/*!
* \ingroup FreeFlowPoreNetworkCoupling
* \brief Coupling conditions for the coupling of a pore-network model
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using FreeFlowPoreNetworkCouplingConditions
    = FreeFlowPoreNetworkCouplingConditionsImplementation<
        MDTraits, CouplingManager,
        GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
        (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)
    >;

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief A base class which provides some common methods used for free-flow/pore-network coupling.
 */
template<class MDTraits, class CouplingManager>
class FreeFlowPoreNetworkCouplingConditionsImplementationBase
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

    using VelocityGradients = StaggeredVelocityGradients;

public:
    static constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    static constexpr auto poreNetworkIndex = CouplingManager::poreNetworkIndex;
private:

    using AdvectionType = GetPropType<SubDomainTypeTag<poreNetworkIndex>, Properties::AdvectionType>;

    static constexpr bool adapterUsed = ModelTraits<poreNetworkIndex>::numFluidPhases() > 1;
    using IndexHelper = FreeFlowPorousMediumCoupling::IndexHelper<freeFlowMassIndex, poreNetworkIndex, FluidSystem<freeFlowMassIndex>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<freeFlowMassIndex>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<poreNetworkIndex>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(FreeFlowPorousMediumCoupling::IsSameFluidSystem<FluidSystem<freeFlowMassIndex>,
                                                                  FluidSystem<poreNetworkIndex>>::value,
                  "All submodels must use the same fluid system");

    using VelocityVector = GlobalPosition<freeFlowMomentumIndex>;

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
     * \brief Returns the momentum flux across the coupling boundary.
     *
     * For the normal momentum coupling, the porous medium side of the coupling condition
     * is evaluated, i.e. -[p n]^pm.
     *
     */
    template<class Context>
    static NumEqVector<freeFlowMomentumIndex> momentumCouplingCondition(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                                        const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                                                        const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars,
                                                                        const Context& context)
    {
        NumEqVector<freeFlowMomentumIndex> momentumFlux(0.0);
        const auto pnmPhaseIdx = couplingPhaseIdx(poreNetworkIndex);
        const auto [pnmPressure, pnmViscosity] = [&]
        {
            for (const auto& scv : scvs(context.fvGeometry))
            {
                if (scv.dofIndex() == context.poreNetworkDofIdx)
                    return std::make_pair(context.elemVolVars[scv].pressure(pnmPhaseIdx), context.elemVolVars[scv].viscosity(pnmPhaseIdx));
            }
            DUNE_THROW(Dune::InvalidStateException, "No coupled scv found");
        }();

        // set p_freeflow = p_PNM
        momentumFlux[scvf.normalAxis()] = pnmPressure;

        // normalize pressure
        momentumFlux[scvf.normalAxis()] -= elemVolVars.gridVolVars().problem().referencePressure(fvGeometry.element(), fvGeometry, scvf);

        // Explicitly account for dv_i/dx_i, which is NOT part of the actual coupling condition. We do it here for convenience so
        // we do not forget to set it in the problem. We assume that the velocity gradient at the boundary towards the interface is the same
        // as the one in the center of the element. TODO check sign
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& frontalInternalScvf = (*scvfs(fvGeometry, scv).begin());
        momentumFlux[scvf.normalAxis()] -= 2*VelocityGradients::velocityGradII(fvGeometry, frontalInternalScvf, elemVolVars) * pnmViscosity;

        // We do NOT consider the inertia term here. If included, Newton convergence decreases drastically and the solution even does not converge to a reference solution.
        // We furthermore assume creeping flow within the boundary layer thus neglecting this term is physically justified.

        momentumFlux[scvf.normalAxis()] *= scvf.directionSign();

        return momentumFlux;
    }

    template<class Context>
    static VelocityVector interfaceThroatVelocity(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                  const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                                  const Context& context)
    {
        const auto& pnmElement = context.fvGeometry.element();
        const auto& pnmFVGeometry = context.fvGeometry;
        const auto& pnmScvf = pnmFVGeometry.scvf(0);
        const auto& pnmElemVolVars = context.elemVolVars;
        const auto& pnmElemFluxVarsCache = context.elemFluxVarsCache;
        const auto& pnmProblem = pnmElemVolVars.gridVolVars().problem();

        const auto pnmPhaseIdx = couplingPhaseIdx(poreNetworkIndex);
        const Scalar area = pnmElemFluxVarsCache[pnmScvf].throatCrossSectionalArea(pnmPhaseIdx);

        // only proceed if area > 0 in order to prevent division by zero (e.g., when the throat was not invaded yet)
        if (area > 0.0)
        {
            using PNMFluxVariables = GetPropType<SubDomainTypeTag<poreNetworkIndex>, Properties::FluxVariables>;
            PNMFluxVariables fluxVars;
            fluxVars.init(pnmProblem, pnmElement, pnmFVGeometry, pnmElemVolVars, pnmScvf, pnmElemFluxVarsCache);

            const Scalar flux = fluxVars.advectiveFlux(pnmPhaseIdx, [pnmPhaseIdx](const auto& volVars){ return volVars.mobility(pnmPhaseIdx);});

            // account for the orientation of the bulk face.
            VelocityVector velocity = (pnmElement.geometry().corner(1) - pnmElement.geometry().corner(0));
            velocity /= velocity.two_norm();
            velocity *= flux / area;

            // TODO: Multiple throats connected to the same pore?
            return velocity;
        }
        else
            return VelocityVector(0.0);
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
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief Coupling conditions specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class FreeFlowPoreNetworkCouplingConditionsImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public FreeFlowPoreNetworkCouplingConditionsImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = FreeFlowPoreNetworkCouplingConditionsImplementationBase<MDTraits, CouplingManager>;
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

    static_assert(GetPropType<SubDomainTypeTag<ParentType::poreNetworkIndex>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<ParentType::poreNetworkIndex>, Properties::ModelTraits>::numFluidPhases(),
                  "Pore-network model must not be compositional");

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the pore-network domain.
     */
    template<class CouplingContext>
    static Scalar massCouplingCondition(Dune::index_constant<ParentType::poreNetworkIndex> domainI,
                                        Dune::index_constant<ParentType::freeFlowMassIndex> domainJ,
                                        const FVElementGeometry<ParentType::poreNetworkIndex>& fvGeometry,
                                        const SubControlVolume<ParentType::poreNetworkIndex>& scv,
                                        const ElementVolumeVariables<ParentType::poreNetworkIndex>& insideVolVars,
                                        const CouplingContext& context)
    {
        Scalar massFlux(0.0);

        for (const auto& c : context)
        {
            // positive values indicate flux into pore-network region
            const Scalar normalFFVelocity = c.velocity * c.scvf.unitOuterNormal();
            const bool pnmIsUpstream = std::signbit(normalFFVelocity);

            const Scalar pnmDensity = insideVolVars[scv].density(couplingPhaseIdx(ParentType::poreNetworkIndex));
            const Scalar ffDensity = c.volVars.density(couplingPhaseIdx(ParentType::freeFlowMassIndex));
            const Scalar area = c.scvf.area() * c.volVars.extrusionFactor();

            // flux is used as source term: positive values mean influx
            massFlux += ParentType::advectiveFlux(pnmDensity, ffDensity, normalFFVelocity*area, pnmIsUpstream);
        }

        return massFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class CouplingContext>
    static Scalar massCouplingCondition(Dune::index_constant<ParentType::freeFlowMassIndex> domainI,
                                        Dune::index_constant<ParentType::poreNetworkIndex> domainJ,
                                        const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                        const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                        const ElementVolumeVariables<ParentType::freeFlowMassIndex>& insideVolVars,
                                        const CouplingContext& context)
    {
        // positive values indicate flux into pore-network region
        const Scalar normalFFVelocity = context.velocity * scvf.unitOuterNormal();
        const bool ffIsUpstream = !std::signbit(normalFFVelocity);

        const Scalar ffDensity = insideVolVars[scvf.insideScvIdx()].density(couplingPhaseIdx(ParentType::freeFlowMassIndex));
        const Scalar pnmDensity = context.volVars.density(couplingPhaseIdx(ParentType::poreNetworkIndex));

        // flux is used in Neumann condition: positive values mean flux out of free-flow domain
        return ParentType::advectiveFlux(ffDensity, pnmDensity, normalFFVelocity, ffIsUpstream);
    }

private:

};

} // end namespace Dumux

#endif
