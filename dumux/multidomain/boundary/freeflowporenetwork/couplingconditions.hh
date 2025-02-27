// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <dune/common/fvector.hh>

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
        GetPropType<typename MDTraits::template SubDomain<1>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
        (GetPropType<typename MDTraits::template SubDomain<1>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)
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

    /*!
     * \brief Evaluate the conductive energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar conductiveEnergyFlux(Dune::index_constant<i> domainI,
                                       Dune::index_constant<j> domainJ,
                                       const SubControlVolumeFace<freeFlowMassIndex>& scvf,
                                       const SubControlVolume<i>& scvI,
                                       const SubControlVolume<j>& scvJ,
                                       const VolumeVariables<i>& volVarsI,
                                       const VolumeVariables<j>& volVarsJ)
    {
        // use properties (distance and thermal conductivity) for transimissibillity coefficient
        //  only from FF side as vertex (pore body) of PNM grid lies on boundary
        const auto& freeFlowVolVars = std::get<const VolumeVariables<freeFlowMassIndex>&>(std::forward_as_tuple(volVarsI, volVarsJ));
        const auto& ffScv = std::get<const SubControlVolume<freeFlowMassIndex>&>(std::forward_as_tuple(scvI, scvJ));
        // distance from FF cell center to interface
        const Scalar distance = getDistance_(ffScv, scvf);
        const Scalar tij = freeFlowVolVars.fluidThermalConductivity() / distance;

        const Scalar deltaT = volVarsJ.temperature() - volVarsI.temperature();

        return -deltaT * tij;
    }

protected:
    /*!
     * \brief Returns the distance between an scvf and the corresponding scv center.
     */
    template<class Scv, class Scvf>
    static Scalar getDistance_(const Scv& scv, const Scvf& scvf)
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
        const auto& pnmVolVars = insideVolVars[scv.indexInElement()];

        for (const auto& c : context)
        {
            // positive values of normal free-flow velocity indicate flux leaving the free flow into the pore-network region
            const Scalar normalFFVelocity = c.velocity * c.scvf.unitOuterNormal();

            // normal pnm velocity (correspondign to its normal vector) is in the opposite direction of normal free flow velocity
            // positive values of normal pnm velocity indicate flux leaving the pnm into the free-flow region
            const Scalar normalPNMVelocity = -normalFFVelocity;
            const bool pnmIsUpstream = std::signbit(normalFFVelocity);
            const Scalar area = c.scvf.area() * c.volVars.extrusionFactor();

            auto flux = massFlux_(domainI, domainJ, c.scvf, scv, c.scv, pnmVolVars, c.volVars, normalPNMVelocity, pnmIsUpstream);
            // flux is used as source term (volumetric flux): positive values mean influx
            // thus, it is multiplied with area and we flip the sign
            flux *= area;
            flux *= -1.0;

            massFlux += flux;
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
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& ffVolVars = insideVolVars[scvf.insideScvIdx()];

        // flux is used in Neumann condition: positive values mean flux out of free-flow domain
        return massFlux_(domainI, domainJ, scvf, insideScv, context.scv, ffVolVars, context.volVars, normalFFVelocity, ffIsUpstream);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the pore network.
     */
    template<class CouplingContext, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyCouplingCondition(Dune::index_constant<ParentType::poreNetworkIndex> domainI,
                                          Dune::index_constant<ParentType::freeFlowMassIndex> domainJ,
                                          const FVElementGeometry<ParentType::poreNetworkIndex>& fvGeometry,
                                          const SubControlVolume<ParentType::poreNetworkIndex>& scv,
                                          const ElementVolumeVariables<ParentType::poreNetworkIndex>& insideVolVars,
                                          const CouplingContext& context)
    {
        Scalar energyFlux(0.0);

        //use VolumeVariables (same type as for context.volVars) instead of ElementVolumeVariables
        const auto& pnmVolVars = insideVolVars[scv.indexInElement()];

        for(const auto& c : context)
        {
            // positive values indicate flux into pore-network region
            const Scalar normalFFVelocity = c.velocity * c.scvf.unitOuterNormal();
            const bool pnmIsUpstream = std::signbit(normalFFVelocity);
            const Scalar normalPNMVelocity = -normalFFVelocity;
            const Scalar area = c.scvf.area() * c.volVars.extrusionFactor();

            auto flux = energyFlux_(domainI, domainJ, c.scvf, scv, c.scv, pnmVolVars, c.volVars, normalPNMVelocity, pnmIsUpstream);

            // flux is used as source term (volumetric flux): positive values mean influx
            // thus, it is multiplied with area and we flip the sign
            flux *= area;
            flux *= -1.0;

            energyFlux += flux;
        }

        return energyFlux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class CouplingContext, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyCouplingCondition(Dune::index_constant<ParentType::freeFlowMassIndex> domainI,
                                          Dune::index_constant<ParentType::poreNetworkIndex> domainJ,
                                          const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                          const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                          const ElementVolumeVariables<ParentType::freeFlowMassIndex>& insideVolVars,
                                          const CouplingContext& context)
    {
        // positive values indicate flux into pore-network region
        const Scalar normalFFVelocity = context.velocity * scvf.unitOuterNormal();
        const bool ffIsUpstream = !std::signbit(normalFFVelocity);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        //use VolumeVariables (same type as for context.volVars) instead of ElementVolumeVariables
        const auto& ffVolVars = insideVolVars[scvf.insideScvIdx()];

        return energyFlux_(domainI, domainJ, scvf, insideScv, context.scv, ffVolVars, context.volVars, normalFFVelocity, ffIsUpstream);
    }

private:
    /*!
     * \brief Evaluate the mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    static Scalar massFlux_(Dune::index_constant<i> domainI,
                                   Dune::index_constant<j> domainJ,
                                   const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                   const SubControlVolume<i>& scvI,
                                   const SubControlVolume<j>& scvJ,
                                   const VolumeVariables<i>& insideVolVars,
                                   const VolumeVariables<j>& outsideVolVars,
                                   const Scalar velocity,
                                   const bool insideIsUpstream)
    {
        Scalar flux(0.0);

        const Scalar insideDensity = insideVolVars.density(couplingPhaseIdx(domainI));
        const Scalar outsideDensity = outsideVolVars.density(couplingPhaseIdx(domainJ));

        flux = ParentType::advectiveFlux(insideDensity, outsideDensity, velocity, insideIsUpstream);

        return flux;
    }

    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyFlux_(Dune::index_constant<i> domainI,
                              Dune::index_constant<j> domainJ,
                              const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                              const SubControlVolume<i>& scvI,
                              const SubControlVolume<j>& scvJ,
                              const VolumeVariables<i>& insideVolVars,
                              const VolumeVariables<j>& outsideVolVars,
                              const Scalar velocity,
                              const bool insideIsUpstream)
    {
        Scalar flux(0.0);

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        flux += ParentType::advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        // conductive energy fluxes
        flux += ParentType::conductiveEnergyFlux(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);

        return flux;
    }
};

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief Coupling conditions specialization for compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class FreeFlowPoreNetworkCouplingConditionsImplementation<MDTraits, CouplingManager, enableEnergyBalance, true>
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
    template<std::size_t id> using FluidSystem = typename VolumeVariables<id>::FluidSystem;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

    static constexpr bool useMoles = GetPropType<SubDomainTypeTag<ParentType::freeFlowMassIndex>, Properties::ModelTraits>::useMoles();
    static constexpr auto numComponents = GetPropType<SubDomainTypeTag<ParentType::freeFlowMassIndex>, Properties::ModelTraits>::numFluidComponents();
    static constexpr auto referenceSystemFormulation = GetPropType<SubDomainTypeTag<ParentType::freeFlowMassIndex>, Properties::MolecularDiffusionType>::referenceSystemFormulation();
    static constexpr auto replaceCompEqIdx = GetPropType<SubDomainTypeTag<ParentType::freeFlowMassIndex>, Properties::ModelTraits>::replaceCompEqIdx();

    static_assert(GetPropType<SubDomainTypeTag<ParentType::poreNetworkIndex>, Properties::ModelTraits>::numFluidComponents() == numComponents,
                  "Models must have same number of components");

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
    using ParentType::couplingCompIdx;

    using NumCompVector = Dune::FieldVector<Scalar, numComponents>;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the pore-network domain.
     */
    template<class CouplingContext>
    static NumCompVector massCouplingCondition(Dune::index_constant<ParentType::poreNetworkIndex> domainI,
                                               Dune::index_constant<ParentType::freeFlowMassIndex> domainJ,
                                               const FVElementGeometry<ParentType::poreNetworkIndex>& fvGeometry,
                                               const SubControlVolume<ParentType::poreNetworkIndex>& scv,
                                               const ElementVolumeVariables<ParentType::poreNetworkIndex>& insideVolVars,
                                               const CouplingContext& context)
    {
        NumCompVector massFlux(0.0);
        const auto& pnmVolVars = insideVolVars[scv.indexInElement()];

        for (const auto& c : context)
        {
            // positive values indicate flux into pore-network region
            const Scalar normalFFVelocity = c.velocity * c.scvf.unitOuterNormal();
            const bool pnmIsUpstream = std::signbit(normalFFVelocity);
            const Scalar normalPNMVelocity = -normalFFVelocity;
            const Scalar area = c.scvf.area() * c.volVars.extrusionFactor();

            auto flux = massFlux_(domainI, domainJ, c.scvf, scv, c.scv, pnmVolVars, c.volVars, normalPNMVelocity, pnmIsUpstream);

            // flux is used as source term (volumetric flux): positive values mean influx
            // thus, it is multiplied with area and we flip the sign
            flux *= area;
            flux *= -1.0;

            massFlux += flux;
        }

        return massFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class CouplingContext>
    static NumCompVector massCouplingCondition(Dune::index_constant<ParentType::freeFlowMassIndex> domainI,
                                               Dune::index_constant<ParentType::poreNetworkIndex> domainJ,
                                               const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                               const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                               const ElementVolumeVariables<ParentType::freeFlowMassIndex>& insideVolVars,
                                               const CouplingContext& context)
    {
        // positive values indicate flux into pore-network region
        const Scalar normalFFVelocity = context.velocity * scvf.unitOuterNormal();
        const bool ffIsUpstream = !std::signbit(normalFFVelocity);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& ffVolVars = insideVolVars[scvf.insideScvIdx()];

        return massFlux_(domainI, domainJ, scvf, insideScv, context.scv, ffVolVars, context.volVars, normalFFVelocity, ffIsUpstream);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the pore network.
     */
    template<class CouplingContext, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyCouplingCondition(Dune::index_constant<ParentType::poreNetworkIndex> domainI,
                                          Dune::index_constant<ParentType::freeFlowMassIndex> domainJ,
                                          const FVElementGeometry<ParentType::poreNetworkIndex>& fvGeometry,
                                          const SubControlVolume<ParentType::poreNetworkIndex>& scv,
                                          const ElementVolumeVariables<ParentType::poreNetworkIndex>& insideVolVars,
                                          const CouplingContext& context)
    {
        Scalar energyFlux(0.0);

        //use VolumeVariables (same type as for context.volVars) instead of ElementVolumeVariables
        const auto& pnmVolVars = insideVolVars[scv.indexInElement()];

        for(const auto& c : context)
        {
            // positive values indicate flux into pore-network region
            const Scalar normalFFVelocity = c.velocity * c.scvf.unitOuterNormal();
            const bool pnmIsUpstream = std::signbit(normalFFVelocity);
            const Scalar normalPNMVelocity = -normalFFVelocity;
            const Scalar area = c.scvf.area() * c.volVars.extrusionFactor();

            auto flux = energyFlux_(domainI, domainJ, c.scvf, scv, c.scv, pnmVolVars, c.volVars, normalPNMVelocity, pnmIsUpstream);

            // flux is used as source term (volumetric flux): positive values mean influx
            // thus, it is multiplied with area and we flip the sign
            flux *= area;
            flux *= -1.0;
            energyFlux += flux;
        }

        return energyFlux;
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class CouplingContext, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyCouplingCondition(Dune::index_constant<ParentType::freeFlowMassIndex> domainI,
                                          Dune::index_constant<ParentType::poreNetworkIndex> domainJ,
                                          const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                          const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                          const ElementVolumeVariables<ParentType::freeFlowMassIndex>& insideVolVars,
                                          const CouplingContext& context)
    {
        // positive values indicate flux into pore-network region
        const Scalar normalFFVelocity = context.velocity * scvf.unitOuterNormal();
        const bool ffIsUpstream = !std::signbit(normalFFVelocity);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        //use VolumeVariables (same type as for context.volVars) instead of ElementVolumeVariables
        const auto& ffVolVars = insideVolVars[scvf.insideScvIdx()];

        return energyFlux_(domainI, domainJ, scvf, insideScv, context.scv, ffVolVars, context.volVars, normalFFVelocity, ffIsUpstream);
    }

private:
    /*!
     * \brief Evaluate the compositional mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    static NumCompVector massFlux_(Dune::index_constant<i> domainI,
                                   Dune::index_constant<j> domainJ,
                                   const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                   const SubControlVolume<i>& scvI,
                                   const SubControlVolume<j>& scvJ,
                                   const VolumeVariables<i>& insideVolVars,
                                   const VolumeVariables<j>& outsideVolVars,
                                   const Scalar velocity,
                                   const bool insideIsUpstream)
    {
        NumCompVector flux(0.0);

        auto moleOrMassFraction = [&](const auto& volVars, int phaseIdx, int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto moleOrMassDensity = [&](const auto& volVars, int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // advective fluxes
        auto insideTerm = [&](int compIdx)
        { return moleOrMassFraction(insideVolVars, couplingPhaseIdx(domainI), compIdx) * moleOrMassDensity(insideVolVars, couplingPhaseIdx(domainI)); };

        auto outsideTerm = [&](int compIdx)
        { return moleOrMassFraction(outsideVolVars, couplingPhaseIdx(domainJ), compIdx) * moleOrMassDensity(outsideVolVars, couplingPhaseIdx(domainJ)); };


        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            flux[domainICompIdx] += ParentType::advectiveFlux(insideTerm(domainICompIdx), outsideTerm(domainJCompIdx), velocity, insideIsUpstream);
        }

        // diffusive fluxes
        NumCompVector diffusiveFlux = diffusiveMolecularFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);

       //convert to correct units if necessary
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged && useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] /= FluidSystem<i>::molarMass(domainICompIdx);
            }
        }
        if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged && !useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] *= FluidSystem<i>::molarMass(domainICompIdx);
            }
        }

        // total flux
        flux += diffusiveFlux;

        // convert to total mass/mole balance, if set be user
        if (replaceCompEqIdx < numComponents)
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);

        return flux;
    }

    /*!
     * \brief Evaluate the energy flux (convective and conductive) across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    static Scalar energyFlux_(Dune::index_constant<i> domainI,
                              Dune::index_constant<j> domainJ,
                              const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                              const SubControlVolume<i>& scvI,
                              const SubControlVolume<j>& scvJ,
                              const VolumeVariables<i>& insideVolVars,
                              const VolumeVariables<j>& outsideVolVars,
                              const Scalar velocity,
                              const bool insideIsUpstream)
    {
        Scalar flux(0.0);

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        flux += ParentType::advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        // conductive fluxes
        flux += ParentType::conductiveEnergyFlux(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);

        // diffusive energy flux
        auto diffusiveFlux = diffusiveMolecularFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);
        for (int compIdx = 0; compIdx < diffusiveFlux.size(); ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            const bool insideDiffFluxIsUpstream = diffusiveFlux[domainICompIdx] > 0.0;
            const Scalar componentEnthalpy = insideDiffFluxIsUpstream ?
                                             getComponentEnthalpy_(insideVolVars, couplingPhaseIdx(domainI), domainICompIdx)
                                           : getComponentEnthalpy_(outsideVolVars, couplingPhaseIdx(domainJ), domainJCompIdx);

           if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
               flux += diffusiveFlux[domainICompIdx] * componentEnthalpy;
           else
               flux += diffusiveFlux[domainICompIdx] * FluidSystem<i>::molarMass(domainICompIdx) * componentEnthalpy;
        }

        return flux;
    }

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    static NumCompVector diffusiveMolecularFlux_(Dune::index_constant<i> domainI,
                                                 Dune::index_constant<j> domainJ,
                                                 const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                                 const SubControlVolume<i>& scvI,
                                                 const SubControlVolume<j>& scvJ,
                                                 const VolumeVariables<i>& volVarsI,
                                                 const VolumeVariables<j>& volVarsJ)
    {
        NumCompVector diffusiveFlux(0.0);
        const Scalar avgDensity = 0.5*(massOrMolarDensity(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI))
                                     + massOrMolarDensity(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ)));

        const auto& freeFlowVolVars = std::get<const VolumeVariables<ParentType::freeFlowMassIndex>&>(std::forward_as_tuple(volVarsI, volVarsJ));
        const auto& ffScv = std::get<const SubControlVolume<ParentType::freeFlowMassIndex>&>(std::forward_as_tuple(scvI, scvJ));
        const Scalar distance = ParentType::getDistance_(ffScv, scvf);

        for (int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            const int freeFlowMainCompIdx = couplingPhaseIdx(ParentType::freeFlowMassIndex);
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            const Scalar massOrMoleFractionI = massOrMoleFraction(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI), domainICompIdx);
            const Scalar massOrMoleFractionJ = massOrMoleFraction(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ), domainJCompIdx);
            const Scalar deltaMassOrMoleFrac = massOrMoleFractionJ - massOrMoleFractionI;

            const Scalar tij = freeFlowVolVars.effectiveDiffusionCoefficient(couplingPhaseIdx(ParentType::freeFlowMassIndex),
                                                                             freeFlowMainCompIdx,
                                                                             couplingCompIdx(ParentType::freeFlowMassIndex, compIdx))
                                                                             / distance;
            diffusiveFlux[domainICompIdx] += -avgDensity * tij * deltaMassOrMoleFrac;
        }

        const Scalar cumulativeFlux = std::accumulate(diffusiveFlux.begin(), diffusiveFlux.end(), 0.0);
        diffusiveFlux[couplingCompIdx(domainI, 0)] = -cumulativeFlux;

        return diffusiveFlux;
    }

    static Scalar getComponentEnthalpy_(const VolumeVariables<ParentType::freeFlowMassIndex>& volVars, int phaseIdx, int compIdx)
    {
        return FluidSystem<ParentType::freeFlowMassIndex>::componentEnthalpy(volVars.fluidState(), 0, compIdx);
    }

    static Scalar getComponentEnthalpy_(const VolumeVariables<ParentType::poreNetworkIndex>& volVars, int phaseIdx, int compIdx)
    {
        return FluidSystem<ParentType::poreNetworkIndex>::componentEnthalpy(volVars.fluidState(), phaseIdx, compIdx);
    }
};

} // end namespace Dumux

#endif
