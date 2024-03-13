// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMConstraint
 */

#ifndef DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_COUPLINGMAPPER_HH
#define DUMUX_MULTIDOMAIN_PORENETWORK_CONSTRAINT_COUPLINGMAPPER_HH

#include <unordered_map>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup PNMConstraint
 * \brief the default mapper for adding a constraint to a throat
 */
template<class MDTraits>
class PNMConstraintCouplingMapper
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
    static constexpr bool isBox()
    { return GridGeometry<i>::discMethod == DiscretizationMethods::box; }

    using MapType = std::unordered_map<std::size_t, std::vector<std::size_t>>;

    static constexpr std::size_t numSD = MDTraits::numSubDomains;

public:
    /*!
     * \brief Main update routine
     */
    template<class CouplingManager>
    void update(const CouplingManager& couplingManager)
    {
        static_assert(numSD == 2, "More than two subdomains not implemented!");
        static_assert(isBox<0>(), "Only box implemented!");

        Dune::Timer watch;
        std::cout << "Initializing the coupling map..." << std::endl;

        for (std::size_t domIdx = 0; domIdx < numSD; ++domIdx)
        {
            stencils_[domIdx].clear();
        }

        const auto& problem0 = couplingManager.problem(domainIdx<0>());
        const auto& problem1 = couplingManager.problem(domainIdx<1>());
        const auto& gg0 = problem0.gridGeometry();
        const auto& gg1 = problem1.gridGeometry();

        for (const auto& element0 : elements(gg0.gridView()))
        {
            auto fvGeometry0 = localView(gg0);
            fvGeometry0.bindElement(element0);

            // add the pair to the map
            const auto eIdx0 = gg0.elementMapper().index(element0);
            // we have the same grid
            const auto eIdx1 = gg1.elementMapper().index(element0);
            stencils_[0][eIdx0].push_back(eIdx1);

            const auto& element1 = gg1.element(eIdx1);
            auto fvGeometry1 = localView(gg1);
            fvGeometry1.bindElement(element1);

            // each element (throat) couples to the connected pore bodies
            for (const auto& scv : scvs(fvGeometry0))
                stencils_[1][eIdx1].push_back(scv.dofIndex());
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
     */
    template<std::size_t i, std::size_t j>
    const std::vector<std::size_t>& couplingStencil(Dune::index_constant<i> domainI,
                                                    const std::size_t eIdxI,
                                                    Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");
        if (isCoupledElement(domainI, eIdxI))
            return stencils_[i].at(eIdxI);

        DUNE_THROW(Dune::InvalidStateException, "Every element needs to be coupled");
    }

    /*!
     * \brief returns true if an element (throat) is coupled
     */
    template<std::size_t i>
    bool isCoupledElement(Dune::index_constant<i>, std::size_t eIdx) const
    { return static_cast<bool>(stencils_[i].count(eIdx)); }

private:
    std::array<MapType, numSD> stencils_;
};

} // end namespace Dumux

#endif
