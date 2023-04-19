// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \copydoc Dumux::FreeFlowPoreNetworkCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORENETWORK_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORENETWORK_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/indices.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/boundary/freeflowporenetwork/ffmassporenetwork/couplingmanager.hh>
#include <dumux/multidomain/boundary/freeflowporenetwork/ffmomentumporenetwork/couplingmanager.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/multibinarycouplingmanager.hh>

#include "couplingconditions.hh"
#include "couplingmapper.hh"

namespace Dumux {

namespace FreeFlowPoreNetworkDetail {

// global subdomain indices
static constexpr auto freeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMassIndex = Dune::index_constant<1>();
static constexpr auto poreNetworkIndex = Dune::index_constant<2>();

// coupling indices
static constexpr auto freeFlowMassToFreeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMomentumToPoreNetworkIndex = Dune::index_constant<1>();
static constexpr auto freeFlowMassToPoreNetworkIndex = Dune::index_constant<2>();
static constexpr auto noCouplingIdx = Dune::index_constant<99>();

constexpr auto makeCouplingManagerMap()
{
    auto map = std::array<std::array<std::size_t, 3>, 3>{};

    // free flow (momentum-mass)
    map[freeFlowMomentumIndex][freeFlowMassIndex] = freeFlowMassToFreeFlowMomentumIndex;
    map[freeFlowMassIndex][freeFlowMomentumIndex] = freeFlowMassToFreeFlowMomentumIndex;

    // free flow momentum - porous medium
    map[freeFlowMomentumIndex][poreNetworkIndex] = freeFlowMomentumToPoreNetworkIndex;
    map[poreNetworkIndex][freeFlowMomentumIndex] = freeFlowMomentumToPoreNetworkIndex;

    // free flow mass - porous medium
    map[freeFlowMassIndex][poreNetworkIndex] = freeFlowMassToPoreNetworkIndex;
    map[poreNetworkIndex][freeFlowMassIndex] = freeFlowMassToPoreNetworkIndex;

    return map;
}

template<std::size_t i>
constexpr auto coupledDomains(Dune::index_constant<i> domainI)
{
    if constexpr (i == freeFlowMomentumIndex)
        return std::make_tuple(freeFlowMassIndex, poreNetworkIndex);
    else if constexpr (i == freeFlowMassIndex)
        return std::make_tuple(freeFlowMomentumIndex, poreNetworkIndex);
    else // i == poreNetworkIndex
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
        return  FreeFlowPoreNetworkDetail::makeCouplingManagerMap();
    }

    template<std::size_t i, std::size_t j>
    static constexpr auto globalToLocal(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        return FreeFlowPoreNetworkDetail::globalToLocalDomainIndices(domainI, domainJ);
    }

    template<std::size_t i>
    static constexpr auto coupledDomains(Dune::index_constant<i> domainI)
    {
        return FreeFlowPoreNetworkDetail::coupledDomains(domainI);
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

    using FreeFlowMomentumPoreNetworkTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMomentumIndex>, SubDomainTypeTag<poreNetworkIndex>
    >;

    using FreeFlowMassPoreNetworkTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMassIndex>, SubDomainTypeTag<poreNetworkIndex>
    >;

    using FreeFlowCouplingManager
        = Dumux::FreeFlowCouplingManager<FreeFlowTraits>;
    using FreeFlowMomentumPoreNetworkCouplingManager
        = Dumux::FreeFlowMomentumPoreNetworkCouplingManager<FreeFlowMomentumPoreNetworkTraits>;
    using FreeFlowMassPoreNetworkCouplingManager
        = Dumux::FreeFlowMassPoreNetworkCouplingManager<FreeFlowMassPoreNetworkTraits>;
};

} // end namespace Detail

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief Coupling manager for coupling freeflow and pore-network models
 */
template<class MDTraits>
class FreeFlowPoreNetworkCouplingManager
: public MultiBinaryCouplingManager<
    MDTraits,
    FreeFlowPoreNetworkDetail::CouplingMaps,
    typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
    typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPoreNetworkCouplingManager,
    typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowMassPoreNetworkCouplingManager
>
{
    using ParentType = MultiBinaryCouplingManager<
        MDTraits,
        FreeFlowPoreNetworkDetail::CouplingMaps,
        typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
        typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPoreNetworkCouplingManager,
        typename FreeFlowPoreNetworkDetail::CouplingManagers<MDTraits>::FreeFlowMassPoreNetworkCouplingManager
    >;

    using ThisType = FreeFlowPoreNetworkCouplingManager<MDTraits>;

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

    using CouplingConditions = FreeFlowPoreNetworkCouplingConditions<MDTraits, FreeFlowPoreNetworkCouplingManager<MDTraits>>;
    using CouplingMapper = StaggeredFreeFlowPoreNetworkCouplingMapper;

public:

    template<std::size_t i, std::size_t j>
    using SubCouplingManager = typename ParentType::template SubCouplingManager<i, j>;

    static constexpr auto freeFlowMomentumIndex = FreeFlowPoreNetworkDetail::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = FreeFlowPoreNetworkDetail::freeFlowMassIndex;
    static constexpr auto poreNetworkIndex = FreeFlowPoreNetworkDetail::poreNetworkIndex;

public:
    using ParentType::ParentType;

    template<class GridVarsTuple>
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<poreNetworkIndex>> poreNetworkProblem,
              GridVarsTuple&& gridVarsTuple,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol); // generic coupling manager stores tuple of shared_ptr

        auto couplingMapper = std::make_shared<CouplingMapper>();
        couplingMapper->update(freeFlowMomentumProblem->gridGeometry(),
            freeFlowMassProblem->gridGeometry(),
            poreNetworkProblem->gridGeometry()
        );

        // initialize the binary sub coupling managers
        typename SubCouplingManager<freeFlowMomentumIndex, freeFlowMassIndex>::SolutionVectorStorage ffSolVecTuple;
        std::get<0>(ffSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).init(
            freeFlowMomentumProblem, freeFlowMassProblem,
            std::make_tuple(std::get<freeFlowMomentumIndex>(gridVarsTuple), std::get<freeFlowMassIndex>(gridVarsTuple)),
            ffSolVecTuple
        );

        typename SubCouplingManager<freeFlowMassIndex, poreNetworkIndex>::SolutionVectorStorage ffMassPmSolVecTuple;
        std::get<0>(ffMassPmSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        std::get<1>(ffMassPmSolVecTuple) = std::get<poreNetworkIndex>(this->curSol());
        this->subCouplingManager(freeFlowMassIndex, poreNetworkIndex).init(
            freeFlowMassProblem, poreNetworkProblem, couplingMapper, ffMassPmSolVecTuple
        );

        typename SubCouplingManager<freeFlowMomentumIndex, poreNetworkIndex>::SolutionVectorStorage ffMomentumPmSolVecTuple;
        std::get<0>(ffMomentumPmSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffMomentumPmSolVecTuple) = std::get<poreNetworkIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, poreNetworkIndex).init(
            freeFlowMomentumProblem, poreNetworkProblem,
            std::make_tuple(std::get<freeFlowMomentumIndex>(gridVarsTuple), std::get<poreNetworkIndex>(gridVarsTuple)),
            couplingMapper, ffMomentumPmSolVecTuple
        );
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary.
     */
    auto massCouplingCondition(Dune::index_constant<poreNetworkIndex> domainI, Dune::index_constant<freeFlowMassIndex> domainJ,
                               const FVElementGeometry<poreNetworkIndex>& fvGeometry,
                               const typename FVElementGeometry<poreNetworkIndex>::SubControlVolume& scv,
                               const ElementVolumeVariables<poreNetworkIndex>& elemVolVars) const
    {

        const auto& couplingContext = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingContext(ii, fvGeometry, scv);
        });

        const auto& freeFlowMassGridGeometry = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.problem(jj).gridGeometry();
        });

        for (auto& c : couplingContext)
        {
            const auto& freeFlowElement = freeFlowMassGridGeometry.element(c.scv.elementIndex());
            c.velocity = this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(freeFlowElement, c.scvf);
        }

        return CouplingConditions::massCouplingCondition(domainI, domainJ, fvGeometry, scv, elemVolVars, couplingContext);
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary.
     */
    auto massCouplingCondition(Dune::index_constant<freeFlowMassIndex> domainI, Dune::index_constant<poreNetworkIndex> domainJ,
                               const FVElementGeometry<freeFlowMassIndex>& fvGeometry,
                               const typename FVElementGeometry<freeFlowMassIndex>::SubControlVolumeFace& scvf,
                               const ElementVolumeVariables<freeFlowMassIndex>& elemVolVars) const
    {

        const auto& couplingContext = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingContext(ii, fvGeometry, scvf);
        });

        couplingContext.velocity = this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(fvGeometry.element(), scvf);
        return CouplingConditions::massCouplingCondition(domainI, domainJ, fvGeometry, scvf, elemVolVars, couplingContext);
    }

    //////////////////////// Conditions for FreeFlowMomentum - PoreNetwork coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    NumEqVector<freeFlowMomentumIndex> momentumCouplingCondition(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                                                 Dune::index_constant<poreNetworkIndex> domainJ,
                                                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                                 const typename FVElementGeometry<freeFlowMomentumIndex>::SubControlVolumeFace& scvf,
                                                                 const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars) const
    {
        if (scvf.isLateral())
            return NumEqVector<freeFlowMomentumIndex>(0.0);

        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, poreNetworkIndex).couplingContext(
            fvGeometry, scvf
        );

        return CouplingConditions::momentumCouplingCondition(fvGeometry, scvf, elemVolVars, context);
    }

    Scalar coupledPoreInscribedRadius(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                      const typename FVElementGeometry<freeFlowMomentumIndex>::SubControlVolumeFace& scvf) const
    {
        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, poreNetworkIndex).couplingContext(
            fvGeometry, scvf
        );

        const auto& pnmScv = [&]
        {
            for (const auto& scv : scvs(context.fvGeometry))
                if (scv.dofIndex() == context.poreNetworkDofIdx)
                    return scv;

            DUNE_THROW(Dune::InvalidStateException, "No scv found");
        }();

        return context.elemVolVars[pnmScv].poreInscribedRadius();
    }

    auto interfaceThroatVelocity(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                 const typename FVElementGeometry<freeFlowMomentumIndex>::SubControlVolumeFace& scvf) const
    {
        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, poreNetworkIndex).couplingContext(
            fvGeometry, scvf
        );

        return CouplingConditions::interfaceThroatVelocity(fvGeometry, scvf, context);
    }

    // //////////////////////// Conditions for FreeFlowMomentum - FreeFlowMass coupling //////////
    // ///////////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).pressure(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the pressure at the center of a sub control volume corresponding to a given sub control volume face.
     *        This is used for setting a Dirichlet pressure for the mass model when a fixed pressure for the momentum balance is set at another
     *        boundary. Since the the pressure at the given scvf is solution-dependent and thus unknown a priori, we just use the value
     *        of the interior cell here.
     */
    Scalar cellPressure(const Element<freeFlowMomentumIndex>& element,
                        const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).cellPressure(
            element, scvf
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    auto insideAndOutsideDensity(const Element<freeFlowMomentumIndex>& element,
                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                 const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                 const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).insideAndOutsideDensity(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, scv, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).effectiveViscosity(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<freeFlowMassIndex>& element,
                      const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(
            element, scvf
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

    // /*!
    //  * \brief Returns whether a given scvf is coupled to the other domain
    //  */
    // bool isCoupledLateralScvf(Dune::index_constant<freeFlowMomentumIndex> domainI,
    //                           Dune::index_constant<poreNetworkIndex> domainJ,
    //                           const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    // {
    //     return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
    //         return cm.isCoupledLateralScvf(ii, scvf);
    //     });
    // }


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
