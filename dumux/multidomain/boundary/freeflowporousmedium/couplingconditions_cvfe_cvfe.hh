// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FFPMCouplingConditionsCvfe
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_CVFE_CVFE_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_CVFE_CVFE_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/traits.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/indexhelper.hh>

namespace Dumux {

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class FFPMCouplingConditionsCvfeImpl;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling conditions for a Darcy model and a (Navier-)Stokes model,
 *        both discretized with control-volume finite element schemes.
 *
 * The porous-medium quantities at the interface are interpolated with the
 * shape functions of the coupled porous-medium element.
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
    template<std::size_t id> using GridDiscretization = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using ElementDiscretization = typename GridDiscretization<id>::LocalView;
    template<std::size_t id> using FluidSystem = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id> using ModelTraits = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;

public:
    static constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;
private:

    static constexpr bool adapterUsed = ModelTraits<porousMediumIndex>::numFluidPhases() > 1;
    using IndexHelper = FreeFlowPorousMediumCoupling::IndexHelper<freeFlowMassIndex, porousMediumIndex, FluidSystem<freeFlowMassIndex>, adapterUsed>;

    static constexpr int enableEnergyBalance = ModelTraits<freeFlowMassIndex>::enableEnergyBalance();
    static_assert(ModelTraits<porousMediumIndex>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(FreeFlowPorousMediumCoupling::IsSameFluidSystem<FluidSystem<freeFlowMassIndex>,
                  FluidSystem<porousMediumIndex>>::value,
                  "All submodels must use the same fluid system");

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
     * \note The permeability is taken from the sub-control volume of the boundary
     *       sub-control-volume face closest to the interpolation point, because the
     *       permeability is defined per sub-control volume.
     */
    template<class IpData, class Context>
    static auto darcyPermeability(const ElementDiscretization<freeFlowMomentumIndex>& elemDisc,
                                  const IpData& ipData,
                                  const Context& context)
    {
        const auto& pmElemDisc = context.elemDisc;
        const auto& pos = ipData.global();

        Scalar minDistance = 0.0;
        std::size_t closestScvfIdx = 0;
        bool foundClosestScvf = false;

        for (const auto& pmScvf : scvfs(pmElemDisc))
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
            DUNE_THROW(Dune::InvalidStateException, "Expected at least one SCVF in porous-medium element discretization");

        return context.vars(pmElemDisc.scv(pmElemDisc.scvf(closestScvfIdx).insideScvIdx())).permeability();
    }

    /*!
     * \brief Returns the porous medium density at a given interpolation point
     */
    template<class ElementDiscretization, class ContextOrElementVars, class IpData>
    static Scalar density(const ElementDiscretization& elemDisc,
                          const ContextOrElementVars& vars,
                          const IpData& ipData)
    {
        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);

        const auto& localBasis = elemDisc.feLocalBasis();

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        localBasis.evaluateFunction(ipData.local(), shapeValues);

        // interpolate density
        Scalar density(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            density += vars[localDof].density(pmPhaseIdx)*shapeValues[localDof.index()][0];

        return density;
    }

    /*!
     * \brief Returns the porous medium pressure at a given interpolation point
     */
    template<class ElementDiscretization, class ContextOrElementVars, class IpData>
    static Scalar pressure(const ElementDiscretization& elemDisc,
                           const ContextOrElementVars& vars,
                           const IpData& ipData)
    {
        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);

        const auto& localBasis = elemDisc.feLocalBasis();

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        localBasis.evaluateFunction(ipData.local(), shapeValues);

        // interpolate pressure
        Scalar pressure(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            pressure += vars[localDof].pressure(pmPhaseIdx)*shapeValues[localDof.index()][0];

        return pressure;
    }

    /*!
     * \brief Evaluate an advective flux across the interface with full upwinding.
     */
    static Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream)
    { return (insideIsUpstream ? insideQuantity : outsideQuantity) * volumeFlow; }
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

    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    static_assert(GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");

    static_assert(!enableEnergyBalance, "Non-isothermal coupling is not implemented for control-volume finite element schemes");

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
};

} // end namespace Dumux

#endif
