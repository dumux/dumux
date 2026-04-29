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

#ifndef DUMUX_MULTIDOMAIN_COUPLINGMAPPER_INTERSECTIONS_HH
#define DUMUX_MULTIDOMAIN_COUPLINGMAPPER_INTERSECTIONS_HH

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
 *        using the grid intersections for coupling. This currently only works for matching grids.
 */
template<class MDTraits, class CouplingManager>
class MapperCoupledMatchingIntersections
{
    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    struct IntersectionInfo
    {
        std::size_t eIdxOutside;
        std::size_t iIdxOutside;
    };

    using FlipIntersectionMapType = std::map<std::pair<std::size_t,std::size_t>, IntersectionInfo>;
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

            for (const auto& intersection1 : intersections(gg1.gridView(), element1))
            {
                // skip all non-boundaries
                if (!intersection1.boundary())
                    continue;

                // get elements intersecting with the intersection center
                // for robustness add epsilon in unit outer normal direction
                const auto& is1Geometry = intersection1.geometry();
                const auto eps = (is1Geometry.center() - is1Geometry.corner(0)).two_norm()*1e-8;
                auto globalPos = is1Geometry.center(); globalPos.axpy(eps, intersection1.centerUnitOuterNormal());
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

                for (const auto& intersection0 : intersections(gg0.gridView(), element0))
                {
                    if (!intersection0.boundary())
                        continue;

                    // TODO: generalize to non-matching grids
                    const auto& is0Geometry = intersection0.geometry();
                    const auto dist = (is1Geometry.center() - is0Geometry.center()).two_norm();
                    if (dist > eps)
                        continue;

                    for(const auto& localDof : localDofs(fvGeometry0, intersection0))
                        isCoupledDof_[0][localDof.dofIndex()] = true;

                    for(const auto& localDof : localDofs(fvGeometry1, intersection1))
                        isCoupledDof_[1][localDof.dofIndex()] = true;

                    std::get<0>(intersectionInfo_)[std::make_pair(eIdx0, static_cast<std::size_t>(intersection0.indexInInside()))]
                        = IntersectionInfo{eIdx1, static_cast<std::size_t>(intersection1.indexInInside())};
                    std::get<1>(intersectionInfo_)[std::make_pair(eIdx1, static_cast<std::size_t>(intersection1.indexInInside()))]
                        = IntersectionInfo{eIdx0, static_cast<std::size_t>(intersection0.indexInInside())};
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
