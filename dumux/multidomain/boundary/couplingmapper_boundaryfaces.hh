// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowMomentumPorousMediumCouplingMapper
 */

#ifndef DUMUX_MULTIDOMAIN_COUPLINGMAPPER_BOUNDARYFACES_HH
#define DUMUX_MULTIDOMAIN_COUPLINGMAPPER_BOUNDARYFACES_HH

#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>

#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief the default mapper for conforming equal dimension boundary coupling between two domains
 *        using the boundary faces for coupling. This currently only works for matching grids.
 */
template<class MDTraits, class CouplingManager>
class MapperCoupledMatchingBoundaryFaces
{
    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    struct FaceInfo
    {
        std::size_t eIdxOutside;
        std::size_t iIdxOutside;
    };

    using FlipIntersectionMapType = std::map<std::pair<std::size_t,std::size_t>, FaceInfo>;
    using MapType = std::map<std::size_t, std::vector<std::size_t>>;

    static constexpr std::size_t numSD = MDTraits::numSubDomains;

public:
    using CouplingStencil = std::vector<std::size_t>;

    /*!
     * \brief Main update routine
     */
    void update(const CouplingManager& couplingManager)
    {
        static_assert(numSD == 2, "More than two subdomains not implemented!");

        Dune::Timer watch;
        std::cout << "Initializing the coupling map..." << std::endl;

        for (std::size_t domIdx = 0; domIdx < numSD; ++domIdx)
            stencils_[domIdx].clear();

        std::get<0>(intersectionInfo_).clear();
        std::get<1>(intersectionInfo_).clear();

        const auto& problem0 = couplingManager.problem(domainIdx<0>());
        const auto& problem1 = couplingManager.problem(domainIdx<1>());
        const auto& gg0 = problem0.gridGeometry();
        const auto& gg1 = problem1.gridGeometry();

        isCoupledDof_[0].resize(gg0.numDofs(), false);
        isCoupledDof_[1].resize(gg1.numDofs(), false);

        auto fvGeometry0 = localView(gg0);
        auto fvGeometry1 = localView(gg1);

        for (const auto& element1 : elements(gg1.gridView()))
        {
            fvGeometry1.bindElement(element1);

            bool elementDofsAdded = false;

            for (const auto& boundaryFace1 : boundaryFaces(fvGeometry1))
            {
                // get elements intersecting with the face center
                // for robustness add epsilon in normal direction
                const auto eps = boundaryFace1.area()*1e-8;
                auto globalPos = boundaryFace1.center(); globalPos.axpy(eps, boundaryFace1.unitOuterNormal());
                const auto indices = intersectingEntities(globalPos, gg0.boundingBoxTree());

                // skip if no intersection was found
                if (indices.empty())
                    continue;

                // sanity check
                if (indices.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Are you sure your sub-domain grids are conformingly discretized on the common interface?");

                // add the pair to the multimap
                const auto eIdx1 = gg1.elementMapper().index(element1);
                const auto eIdx0 = indices[0];
                const auto& element0 = gg0.element(eIdx0);
                fvGeometry0.bindElement(element0);

                if(!elementDofsAdded)
                {
                    for(auto&& localDof : localDofs(fvGeometry1))
                        stencils_[0][eIdx0].push_back(localDof.dofIndex());

                    for(auto&& localDof : localDofs(fvGeometry0))
                        stencils_[1][eIdx1].push_back(localDof.dofIndex());

                    elementDofsAdded = true;
                }

                for (const auto& boundaryFace0 : boundaryFaces(fvGeometry0))
                {
                    // TODO: generalize to non-matching grids
                    const auto dist = (boundaryFace1.center() - boundaryFace0.center()).two_norm();
                    if (dist > eps)
                        continue;

                    for(const auto& localDof : localDofs(fvGeometry0, boundaryFace0))
                        isCoupledDof_[0][localDof.dofIndex()] = true;

                    for(const auto& localDof : localDofs(fvGeometry1, boundaryFace1))
                        isCoupledDof_[1][localDof.dofIndex()] = true;

                    std::get<0>(intersectionInfo_)[std::make_pair(eIdx0, static_cast<std::size_t>(boundaryFace0.intersectionIndex()))]
                        = FaceInfo{eIdx1, static_cast<std::size_t>(boundaryFace1.intersectionIndex())};
                    std::get<1>(intersectionInfo_)[std::make_pair(eIdx1, static_cast<std::size_t>(boundaryFace1.intersectionIndex()))]
                        = FaceInfo{eIdx0, static_cast<std::size_t>(boundaryFace0.intersectionIndex())};
                }
            }
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param eIdxI the index of the coupled element of domain i
     * \param domainJ the domain index of domain j
     *
     * \return the coupling stencil, i.e. the indices of the coupled dofs of domain j
     */
    template<std::size_t i, std::size_t j>
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<i> domainI,
                                                    const std::size_t eIdxI,
                                                    Dune::index_constant<j> domainJ) const
    {
        if (isCoupledElement(domainI, eIdxI))
            return stencils_[i].at(eIdxI);
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
     * \brief Return if a dof with index dofIdx of domain i is coupled to domain j
     */
    template<std::size_t i>
    bool isCoupledDof(Dune::index_constant<i>, std::size_t dofIdx) const
    {
        return isCoupledDof_[i][dofIdx];
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param eIdxI the element index
     * \param iIdxI the intersection index
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   std::size_t eIdxI,
                   std::size_t iIdxI) const
    {
        return static_cast<bool>(std::get<i>(intersectionInfo_).count(std::make_pair(eIdxI, iIdxI)));
    }

    /*!
     * \brief Return the outside element index (the element index of the other domain)
     * \param domainI the domain index for which to compute the flux
     * \param eIdxI the element index
     * \param iIdxI the intersection index
     */
    template<std::size_t i>
    std::size_t outsideElementIndex(Dune::index_constant<i> domainI,
                                    std::size_t eIdxI,
                                    std::size_t iIdxI) const
    {
        return std::get<i>(intersectionInfo_).at(std::make_pair(eIdxI, iIdxI)).eIdxOutside;
    }

    /*!
     * \brief Return the intersection index on the outside element (the other domain)
     * \param domainI the domain index for which to compute the flux
     * \param eIdxI the element index
     * \param iIdxI the intersection index
     */
    template<std::size_t i>
    std::size_t flipIntersectionIndex(Dune::index_constant<i> domainI,
                                      std::size_t eIdxI,
                                      std::size_t iIdxI) const
    {
        return std::get<i>(intersectionInfo_).at(std::make_pair(eIdxI, iIdxI)).iIdxOutside;
    }

private:
    std::array<MapType, numSD> stencils_;
    std::vector<std::size_t> emptyStencil_;
    std::array<std::vector<bool>, numSD> isCoupledDof_;
    std::tuple<FlipIntersectionMapType, FlipIntersectionMapType> intersectionInfo_;
};

} // end namespace Dumux

#endif
