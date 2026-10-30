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

#include <cstddef>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/traits.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/indexhelper.hh>

namespace Dumux {

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class FFPMCouplingConditionsCvfeImpl;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling conditions for a Darcy model and a (Navier-)Stokes model,
 *        both discretized with control-volume finite element or finite element schemes.
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
    template<std::size_t id> using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
    template<std::size_t id> using AdvectionType = GetPropType<SubDomainTypeTag<id>, Properties::AdvectionType>;
    template<std::size_t id> using GlobalPosition = typename GridDiscretization<id>::GridView::template Codim<0>::Geometry::GlobalCoordinate;

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
     * \brief Returns the intrinsic permeability of the coupled Darcy element at an interpolation point.
     * \note The permeability is interpolated from the local dofs.
     */
    template<class IpData, class Context>
    static auto darcyPermeability(const ElementDiscretization<freeFlowMomentumIndex>& elemDisc,
                                  const IpData& ipData,
                                  const Context& context)
    {
        return interpolate_(context, ipData.global(), [](const auto& vars) { return vars.permeability(); });
    }

    /*!
     * \brief Returns the Darcy velocity at an interpolation point of the coupled Darcy element.
     * \note The velocity is the flux term of the advection type of the porous-medium model with
     *       the tensor \f$\mathbf{K}/\mu\f$, which is interpolated from the local dofs.
     */
    template<class Context, class IpData>
    static auto darcyVelocity(const Context& context, const IpData& ipData)
    {
        static_assert(!Concept::FVGridVariables<GridVariables<porousMediumIndex>>,
                      "The Darcy velocity requires the porous-medium model to use the general grid variables");

        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);
        const auto& elemDisc = context.elemDisc;

        // the coupling context provides the variables per local dof like the element variables,
        // but holds no interpolation point cache, so the cache is evaluated here
        typename Concept::GridVariablesCache_t<GridVariables<porousMediumIndex>>::InterpolationPointData ipCache;
        ipCache.update(context.problem, elemDisc.element(), elemDisc, context, ipData);

        const auto mobilityTensor = [&](const auto& pmElemDisc, const auto& pmElemVars, const auto& pmIpData)
        {
            return interpolate_(context, pmIpData.global(), [&](const auto& vars)
            {
                auto tensor = vars.permeability();
                tensor /= vars.viscosity(pmPhaseIdx);
                return tensor;
            });
        };

        return AdvectionType<porousMediumIndex>::fluxTerm(context.problem, elemDisc.element(), elemDisc,
                                                          context, pmPhaseIdx, ipData, ipCache, mobilityTensor);
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

private:
    //! a quantity of the porous-medium variables interpolated from the local dofs at a position on the coupling interface
    template<class Context, class Quantity>
    static auto interpolate_(const Context& context,
                             const GlobalPosition<porousMediumIndex>& globalPos,
                             const Quantity& quantity)
    {
        const auto& pmElemDisc = context.elemDisc;

        // the position lies on the face of the free-flow element, so its local position
        // in the coupled porous-medium element is needed
        const auto pmLocalPos = pmElemDisc.elementGeometry().local(globalPos);

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        pmElemDisc.feLocalBasis().evaluateFunction(pmLocalPos, shapeValues);

        std::decay_t<decltype(quantity(context.vars(*localDofs(pmElemDisc).begin())))> value(0.0);
        for (const auto& localDof : localDofs(pmElemDisc))
        {
            auto localValue = quantity(context.vars(localDof));
            localValue *= shapeValues[localDof.index()][0];
            value += localValue;
        }
        return value;
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
