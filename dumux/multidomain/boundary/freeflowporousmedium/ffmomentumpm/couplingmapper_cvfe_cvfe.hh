// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMomentumPorousMediumCouplingMapper
 */

#ifndef DUMUX_MULTIDOMAIN_FREEFLOWMOMENTUM_POROUSMEDIUM_COUPLINGMAPPER_CVFE_CVFE_HH
#define DUMUX_MULTIDOMAIN_FREEFLOWMOMENTUM_POROUSMEDIUM_COUPLINGMAPPER_CVFE_CVFE_HH

#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief the default mapper for conforming equal dimension boundary coupling between two domains (box or cc)
 * \todo how to extend to arbitrary number of domains?
 */
template<class MDTraits, class CouplingManager>
class FFMomentumPMCouplingMapperCvfe
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using SubControlVolume = typename GridGeometry<i>::SubControlVolume;
    template<std::size_t i> using SubControlVolumeFace = typename GridGeometry<i>::SubControlVolumeFace;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;
    template<std::size_t i> using Element = typename GridView<i>::template Codim<0>::Entity;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    template<std::size_t i>
    static constexpr bool isCvfe()
    { return GridGeometry<i>::discMethod == DiscretizationMethods::isCVFE<GridGeometry<i>::discMethod>; }

    struct ScvfInfoPM
    {
        std::size_t eIdxOutside;
        std::size_t flipScvfIdx;
    };

    struct ScvfInfoFF
    {
        std::size_t eIdxOutside;
        std::size_t flipScvfIdx;
    };

    using FlipScvfMapTypePM = std::unordered_map<std::pair<std::size_t,std::size_t>, ScvfInfoPM>;
    using FlipScvfMapTypeFF = std::unordered_map<std::pair<std::size_t,std::size_t>, ScvfInfoFF>;
    using MapType = std::unordered_map<std::size_t, std::vector<std::size_t>>;

    static constexpr std::size_t numSD = MDTraits::numSubDomains;

public:
    /*!
     * \brief Main update routine
     */
    void update(const CouplingManager& couplingManager)
    {
        static_assert(numSD == 2, "More than two subdomains not implemented!");
        static_assert(isCvfe()<0>() && isCvfe()<1>(), "Only coupling between cvfe and cvfe schemes implemented!");

        Dune::Timer watch;
        std::cout << "Initializing the coupling map..." << std::endl;

        for (std::size_t domIdx = 0; domIdx < numSD; ++domIdx)
            stencils_[domIdx].clear();

        std::get<CouplingManager::freeFlowMomentumIndex>(scvfInfo_).clear();
        std::get<CouplingManager::porousMediumIndex>(scvfInfo_).clear();

        const auto& freeFlowMomentumProblem = couplingManager.problem(CouplingManager::freeFlowMomentumIndex);
        const auto& porousMediumProblem = couplingManager.problem(CouplingManager::porousMediumIndex);
        const auto& freeFlowMomentumGG = freeFlowMomentumProblem.gridGeometry();
        const auto& porousMediumGG = porousMediumProblem.gridGeometry();

        isCoupledDof[CouplingManager::freeFlowMomentumIndex]_.resize(freeFlowMomentumGG.numDofs(), false);
        isCoupledDof_[CouplingManager::porousMediumIndex].resize(porousMediumGG.numDofs(), false);

        auto pmFvGeometry = localView(porousMediumGG);
        auto ffFvGeometry = localView(freeFlowMomentumGG);

        for (const auto& pmElement : elements(porousMediumGG.gridView()))
        {
            pmFvGeometry.bindElement(pmElement);

            bool elementDofsAdded = false;

            for (const auto& pmScvf : scvfs(pmFvGeometry))
            {
                // skip all non-boundaries
                if (!pmScvf.boundary())
                    continue;

                // get elements intersecting with the scvf center
                // for robustness add epsilon in unit outer normal direction
                const auto eps = (pmScvf.ipGlobal() - pmFvGeometry.geometry(pmScvf).corner(0)).two_norm()*1e-8;
                auto globalPos = pmScvf.ipGlobal(); globalPos.axpy(eps, pmScvf.unitOuterNormal());
                const auto indices = intersectingEntities(globalPos, freeFlowMomentumGG.boundingBoxTree());

                // skip if no intersection was found
                if (indices.empty())
                    continue;

                // sanity check
                if (indices.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Are you sure your sub-domain grids are conformingly discretized on the common interface?");

                // add the pair to the multimap
                const auto pmElemIdx = porousMediumGG.elementMapper().index(pmElement);
                const auto ffElemIdx = indices[0];
                const auto& ffElement = freeFlowMomentumGG.element(ffElemIdx);
                ffFvGeometry.bindElement(ffElement);

                if(!elementDofsAdded)
                {
                    for(auto&& scv : scvs(pmFvGeometry))
                        stencils_[CouplingManager::freeFlowMomentumIndex][ffElemIdx].push_back(scv.dofIndex());

                    for(auto&& scv : scvs(ffFvGeometry))
                        stencils_[CouplingManager::porousMediumIndex][pmElemIdx].push_back(scv.dofIndex());

                    elementDofsAdded = true;
                }

                for (const auto& ffScvf : scvfs(ffFvGeometry))
                {
                    if (!ffScvf.boundary())
                        continue;

                    // TODO: generalize to non-matching grids
                    const auto dist = (pmScvf.ipGlobal() - ffScvf.ipGlobal()).two_norm();
                    if (dist > eps)
                        continue;

                    isCoupledDof_[CouplingManager::freeFlowMomentumIndex][ffFvGeometry.scv(ffScvf.insideScvIdx()).dofIndex()] = true;
                    isCoupledDof_[CouplingManager::porousMediumIndex][pmFvGeometry.scv(pmScvf.insideScvIdx()).dofIndex()] = true;

                    std::get<CouplingManager::porousMediumIndex>(scvfInfo_)[std::make_pair<std::size_t, std::size_t>(pmElemIdx, pmScvf.index())]
                        = ScvfInfoPM{ffElemIdx, ffScvf.index()};
                    std::get<CouplingManager::freeFlowMomentumIndex>(scvfInfo_)[std::make_pair<std::size_t, std::size_t>(ffElemIdx, ffScvf.index())]
                        = ScvfInfoFF{pmElemIdx, pmScvf.index()};
                }
            }
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param eIdxI the index of the coupled element of domain í
     * \param domainJ the domain index of domain j
     *
     * \note  The element residual definition depends on the discretization scheme of domain i
     *        box: a container of the residuals of all sub control volumes
     *        cc : the residual of the (sub) control volume
     *        fem: the residual of the element
     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<CouplingManager::porousMediumIndex> domainI,
                                                    const std::size_t eIdxI,
                                                    Dune::index_constant<CouplingManager::freeFlowMomentumIndex> domainJ) const
    {
        if (isCoupledElement(domainI, eIdxI))
            return stencils_[CouplingManager::porousMediumIndex].at(eIdxI);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return if an element residual with index eIdx of domain i is coupled to domain j
     */
    template<std::size_t i>
    bool isCoupledElement(Dune::index_constant<i>, std::size_t eIdx) const
    {
        return static_cast<bool>(stencils_[i].count(eIdx));
    }

    /*!
     * \brief Return if an dof with index dofIdx of domain i is coupled to domain j
     */
    template<std::size_t i>
    bool isCoupledDof(Dune::index_constant<i>, std::size_t dofIdx) const
    {
        return isCoupledDof_[i][dofIdx];
    }


    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   std::size_t eIdxI,
                   std::size_t localScvfIdxI) const
    {
        return static_cast<bool>(std::get<i>(scvfInfo_).count(std::make_pair<std::size_t, std::size_t>(eIdxI, localScvfIdxI)));
    }

    /*!
     * \brief Return the scvf index of the flipped scvf in the other domain
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    std::size_t flipScvfIndex(Dune::index_constant<i> domainI,
                              std::size_t eIdxI,
                              std::size_t localScvfIdxI) const
    {
        return std::get<i>(scvfInfo_).at(std::make_pair<std::size_t, std::size_t>(eIdxI, localScvfIdxI)).flipScvfIdx;
    }

    /*!
     * \brief Return the outside element index (the element index of the other domain)
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    std::size_t outsideElementIndex(Dune::index_constant<i> domainI,
                                    std::size_t eIdxI,
                                    std::size_t localScvfIdxI) const
    {
        return std::get<i>(scvfInfo_).at(std::make_pair<std::size_t, std::size_t>(eIdxI, localScvfIdxI)).eIdxOutside;
    }

private:
    std::array<MapType, numSD> stencils_;
    std::vector<std::size_t> emptyStencil_;
    std::array<std::vector<bool>, numSD> isCoupledDof_;
    std::tuple<FlipScvfMapTypeFF, FlipScvfMapTypePM> scvfInfo_;

};

} // end namespace Dumux

#endif
