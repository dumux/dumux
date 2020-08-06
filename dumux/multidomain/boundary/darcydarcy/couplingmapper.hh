// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup DarcyDarcyCoupling
 * \copydoc Dumux::DarcyDarcyBoundaryCouplingMapper
 */

#ifndef DUMUX_MULTIDOMAIN_DARCYDARCY_BOUNDARY_COUPLINGMAPPER_HH
#define DUMUX_MULTIDOMAIN_DARCYDARCY_BOUNDARY_COUPLINGMAPPER_HH

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
 * \ingroup DarcyDarcyCoupling
 * \brief the default mapper for conforming equal dimension boundary coupling between two domains (box or cc)
 * \todo how to extend to arbitrary number of domains?
 */
template<class MDTraits>
class DarcyDarcyBoundaryCouplingMapper
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using SubControlVolumeFace = typename GridGeometry<i>::SubControlVolumeFace;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;
    template<std::size_t i> using Element = typename GridView<i>::template Codim<0>::Entity;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    template<std::size_t i>
    static constexpr bool isCCTpfa()
    { return GridGeometry<i>::discMethod == DiscretizationMethod::cctpfa; }

    struct ScvfInfo
    {
        std::size_t eIdxOutside;
        std::size_t flipScvfIdx;
    };

    using FlipScvfMapType = std::unordered_map<std::size_t, ScvfInfo>;
    using MapType = std::unordered_map<std::size_t, std::vector<std::size_t>>;

    static constexpr std::size_t numSD = MDTraits::numSubDomains;

public:
    /*!
     * \brief Main update routine
     */
    template<class CouplingManager>
    void update(const CouplingManager& couplingManager)
    {
        // TODO: Box and multiple domains
        static_assert(numSD == 2, "More than two subdomains not implemented!");
        static_assert(isCCTpfa<0>() && isCCTpfa<1>(), "Only cctpfa implemented!");

        Dune::Timer watch;
        std::cout << "Initializing the coupling map..." << std::endl;

        for (std::size_t domIdx = 0; domIdx < numSD; ++domIdx)
        {
            stencils_[domIdx].clear();
            scvfInfo_[domIdx].clear();
        }

        const auto& problem0 = couplingManager.problem(domainIdx<0>());
        const auto& problem1 = couplingManager.problem(domainIdx<1>());
        const auto& gg0 = problem0.gridGeometry();
        const auto& gg1 = problem1.gridGeometry();

        isCoupledScvf_[0].resize(gg0.numScvf(), false);
        isCoupledScvf_[1].resize(gg1.numScvf(), false);

        for (const auto& element0 : elements(gg0.gridView()))
        {
            auto fvGeometry0 = localView(gg0);
            fvGeometry0.bindElement(element0);

            for (const auto& scvf0 : scvfs(fvGeometry0))
            {
                // skip all non-boundaries
                if (!scvf0.boundary())
                    continue;

                // get elements intersecting with the scvf center
                // for robustness add epsilon in unit outer normal direction
                const auto eps = (scvf0.ipGlobal() - element0.geometry().corner(0)).two_norm()*1e-8;
                auto globalPos = scvf0.ipGlobal(); globalPos.axpy(eps, scvf0.unitOuterNormal());
                const auto indices = intersectingEntities(globalPos, gg1.boundingBoxTree());

                // skip if no intersection was found
                if (indices.empty())
                    continue;

                // sanity check
                if (indices.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Are you sure your grids is conforming at the boundary?");

                // add the pair to the multimap
                const auto eIdx0 = gg0.elementMapper().index(element0);
                const auto eIdx1 = indices[0];
                stencils_[0][eIdx0].push_back(eIdx1);

                // mark the scvf and find and mark the flip scvf
                isCoupledScvf_[0][scvf0.index()] = true;
                const auto& element1 = gg1.element(eIdx1);
                auto fvGeometry1 = localView(gg1);
                fvGeometry1.bindElement(element1);

                using std::abs;
                for (const auto& scvf1 : scvfs(fvGeometry1))
                    if (scvf1.boundary())
                        if (abs(scvf1.unitOuterNormal()*scvf0.unitOuterNormal() + 1) < eps)
                        {
                            isCoupledScvf_[1][scvf1.index()] = true;
                            scvfInfo_[0][scvf0.index()] = ScvfInfo{eIdx1, scvf1.index()};
                            scvfInfo_[1][scvf1.index()] = ScvfInfo{eIdx0, scvf0.index()};
                        }
            }
        }

        // create the inverse map for efficient access
        for (const auto& entry : stencils_[0])
            for (const auto idx : entry.second)
                stencils_[1][idx].push_back(entry.first);

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
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
    template<std::size_t i, std::size_t j>
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<i> domainI,
                                                    const std::size_t eIdxI,
                                                    Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");
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
    { return static_cast<bool>(stencils_[i].count(eIdx)); }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const SubControlVolumeFace<i>& scvf) const
    {
        return isCoupledScvf_[i].at(scvf.index());
    }

    /*!
     * \brief Return the scvf index of the flipped scvf in the other domain
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    std::size_t flipScvfIndex(Dune::index_constant<i> domainI,
                              const SubControlVolumeFace<i>& scvf) const
    {
        return scvfInfo_[i].at(scvf.index()).flipScvfIdx;
    }

    /*!
     * \brief Return the outside element index (the element index of the other domain)
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    std::size_t outsideElementIndex(Dune::index_constant<i> domainI,
                                    const SubControlVolumeFace<i>& scvf) const
    {
        return scvfInfo_[i].at(scvf.index()).eIdxOutside;
    }

private:
    std::array<MapType, numSD> stencils_;
    std::vector<std::size_t> emptyStencil_;
    std::array<std::vector<bool>, numSD> isCoupledScvf_;
    std::array<FlipScvfMapType, numSD> scvfInfo_;
};

} // end namespace Dumux

#endif
