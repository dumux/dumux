// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Base class for coupling freeflow and porous medium flow models.
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_BASE_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_BASE_HH

#include <utility>
#include <memory>

#include <dune/common/indices.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/ffmasspm/couplingmanager.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/ffmomentumpm/couplingmanager.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/multibinarycouplingmanager.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace FreeFlowPorousMediumDetail {

// global subdomain indices
static constexpr auto freeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMassIndex = Dune::index_constant<1>();
static constexpr auto porousMediumIndex = Dune::index_constant<2>();

// coupling indices
static constexpr auto freeFlowMassToFreeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMomentumToPorousMediumIndex = Dune::index_constant<1>();
static constexpr auto freeFlowMassToPorousMediumIndex = Dune::index_constant<2>();
static constexpr auto noCouplingIdx = Dune::index_constant<99>();

constexpr auto makeCouplingManagerMap()
{
    auto map = std::array<std::array<std::size_t, 3>, 3>{};

    // free flow (momentum-mass)
    map[freeFlowMomentumIndex][freeFlowMassIndex] = freeFlowMassToFreeFlowMomentumIndex;
    map[freeFlowMassIndex][freeFlowMomentumIndex] = freeFlowMassToFreeFlowMomentumIndex;

    // free flow momentum - porous medium
    map[freeFlowMomentumIndex][porousMediumIndex] = freeFlowMomentumToPorousMediumIndex;
    map[porousMediumIndex][freeFlowMomentumIndex] = freeFlowMomentumToPorousMediumIndex;

    // free flow mass - porous medium
    map[freeFlowMassIndex][porousMediumIndex] = freeFlowMassToPorousMediumIndex;
    map[porousMediumIndex][freeFlowMassIndex] = freeFlowMassToPorousMediumIndex;

    return map;
}

template<std::size_t i>
constexpr auto coupledDomains(Dune::index_constant<i> domainI)
{
    if constexpr (i == freeFlowMomentumIndex)
        return std::make_tuple(freeFlowMassIndex, porousMediumIndex);
    else if constexpr (i == freeFlowMassIndex)
        return std::make_tuple(freeFlowMomentumIndex, porousMediumIndex);
    else // i == porousMediumIndex
        return std::make_tuple(freeFlowMomentumIndex, freeFlowMassIndex);
}

template<std::size_t i, std::size_t j>
constexpr auto globalToLocalDomainIndices(Dune::index_constant<i>, Dune::index_constant<j>)
{
    static_assert(i <= 2 && j <= 2);
    static_assert(i != j);

    if constexpr (i < j)
        return std::pair<Dune::index_constant<0>, Dune::index_constant<1>>{};
    else
        return std::pair<Dune::index_constant<1>, Dune::index_constant<0>>{};
}

struct CouplingMaps
{
    static constexpr auto managerMap()
    {
        return  FreeFlowPorousMediumDetail::makeCouplingManagerMap();
    }

    template<std::size_t i, std::size_t j>
    static constexpr auto globalToLocal(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        return FreeFlowPorousMediumDetail::globalToLocalDomainIndices(domainI, domainJ);
    }

    template<std::size_t i>
    static constexpr auto coupledDomains(Dune::index_constant<i> domainI)
    {
        return FreeFlowPorousMediumDetail::coupledDomains(domainI);
    }
};

template<class MDTraits>
struct CouplingManagers
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    using FreeFlowTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMomentumIndex>, SubDomainTypeTag<freeFlowMassIndex>
    >;

    using FreeFlowMomentumPorousMediumTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMomentumIndex>, SubDomainTypeTag<porousMediumIndex>
    >;

    using FreeFlowMassPorousMediumTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMassIndex>, SubDomainTypeTag<porousMediumIndex>
    >;

    using FreeFlowCouplingManager
        = Dumux::FreeFlowCouplingManager<FreeFlowTraits>;
    using FreeFlowMomentumPorousMediumCouplingManager
        = Dumux::FreeFlowMomentumPorousMediumCouplingManager<FreeFlowMomentumPorousMediumTraits>;
    using FreeFlowMassPorousMediumCouplingManager
        = Dumux::FreeFlowMassPorousMediumCouplingManager<FreeFlowMassPorousMediumTraits>;
};

} // end namespace FreeFlowPorousMediumDetail
#endif // DOXYGEN

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Base coupling manager for coupling freeflow and porous medium flow models
 */
template<class MDTraits>
class FreeFlowPorousMediumCouplingManagerBase
: public MultiBinaryCouplingManager<
    MDTraits,
    FreeFlowPorousMediumDetail::CouplingMaps,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPorousMediumCouplingManager,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMassPorousMediumCouplingManager
>
{
    using ParentType = MultiBinaryCouplingManager<
        MDTraits,
        FreeFlowPorousMediumDetail::CouplingMaps,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPorousMediumCouplingManager,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMassPorousMediumCouplingManager
    >;

    using Scalar = typename MDTraits::Scalar;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    using SolutionVector = typename MDTraits::SolutionVector;

public:

    template<std::size_t i, std::size_t j>
    using SubCouplingManager = typename ParentType::template SubCouplingManager<i, j>;

    static constexpr auto freeFlowMomentumIndex = FreeFlowPorousMediumDetail::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = FreeFlowPorousMediumDetail::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = FreeFlowPorousMediumDetail::porousMediumIndex;

public:
    using ParentType::ParentType;

    template<class GridVarsTuple>
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<porousMediumIndex>> porousMediumProblem,
              GridVarsTuple&& gridVarsTuple,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol); // generic coupling manager stores tuple of shared_ptr

        // initialize the binary sub coupling managers
        typename SubCouplingManager<freeFlowMomentumIndex, freeFlowMassIndex>::SolutionVectorStorage ffSolVecTuple;
        std::get<0>(ffSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).init(
            freeFlowMomentumProblem, freeFlowMassProblem,
            std::make_tuple(std::get<freeFlowMomentumIndex>(gridVarsTuple), std::get<freeFlowMassIndex>(gridVarsTuple)),
            ffSolVecTuple
        );

        typename SubCouplingManager<freeFlowMassIndex, porousMediumIndex>::SolutionVectorStorage ffMassPmSolVecTuple;
        std::get<0>(ffMassPmSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        std::get<1>(ffMassPmSolVecTuple) = std::get<porousMediumIndex>(this->curSol());
        this->subCouplingManager(freeFlowMassIndex, porousMediumIndex).init(
            freeFlowMassProblem, porousMediumProblem, ffMassPmSolVecTuple
        );

        typename SubCouplingManager<freeFlowMomentumIndex, porousMediumIndex>::SolutionVectorStorage ffMomentumPmSolVecTuple;
        std::get<0>(ffMomentumPmSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffMomentumPmSolVecTuple) = std::get<porousMediumIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).init(
            freeFlowMomentumProblem, porousMediumProblem, ffMomentumPmSolVecTuple
        );
    }

    template<std::size_t i>
    const Problem<i>& problem(Dune::index_constant<i> domainI) const
    {
        return this->subApply(domainI, [&](const auto& cm, auto&& ii) -> const auto& {
            return cm.problem(ii);
        });
    }

    template<std::size_t i, std::size_t j>
    bool isCoupled(Dune::index_constant<i> domainI,
                   Dune::index_constant<j> domainJ,
                   const SubControlVolumeFace<i>& scvf) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupled(ii, scvf);
        });
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index of domain i for which to compute the flux
     * \param domainJ the domain index of domain j for which to compute the flux
     * \param scv the sub control volume
     */
    template<std::size_t i, std::size_t j>
    bool isCoupled(Dune::index_constant<i> domainI,
                   Dune::index_constant<j> domainJ,
                   const SubControlVolume<i>& scv) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupled(ii, scv);
        });
    }

    using ParentType::couplingStencil;
    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given sub-control volume of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    template<std::size_t j>
    const auto& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                const Element<freeFlowMomentumIndex>& elementI,
                                const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                Dune::index_constant<j> domainJ) const
    {
        static_assert(freeFlowMomentumIndex != j);
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingStencil(ii, elementI, scvI, jj);
        });
    }
};

} // end namespace Dumux

#endif
