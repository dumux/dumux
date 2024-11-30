// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Class to represent a decomposition of a domain into multiple subdomains with mortars on their interfaces.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH
#define DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH

#include <type_traits>
#include <concepts>
#include <variant>
#include <memory>
#include <vector>
#include <ranges>

#include <dune/common/exceptions.hh>
#include <dumux/multidomain/glue.hh>

namespace Dumux::Mortar {

template<typename MortarGridGeometry, typename... GridGeometries>
class DecompositionFactory;

/*!
 * \ingroup MultiDomain
 * \brief Class to represent a decomposition of a domain into multiple subdomains with mortars on their interfaces.
 * \tparam MortarGridGeometry Grid geometry type of the mortar domain.
 * \tparam GridGeometries Types of all subdomain grid geometries that may occur in the decomposition.
 */
template<typename MortarGridGeometry, typename... GridGeometries>
class Decomposition
{
 public:
    //! Visit the mortar domains that are coupled to the given subdomain
    template<typename GridGeometry, std::invocable<const std::shared_ptr<const MortarGridGeometry>&> Visitor>
        requires(std::disjunction_v<std::is_same<GridGeometry, GridGeometries>...>)
    void visitCoupledMortarsOf(const GridGeometry& subDomain, Visitor&& v) const
    {
        std::ranges::for_each(
            subDomainsToMortar_.at(subDomainIndexOf_(subDomain)),
            [&] (auto mortarIdx) { v(mortars_.at(mortarIdx)); }
        );
    }

    //! Visit the subdomains that are coupled to the given mortar domain
    template<typename Visitor>
        requires(std::conjunction_v<std::is_invocable<Visitor, const std::shared_ptr<const GridGeometries>&>...>)
    void visitCoupledSubDomainsOf(const MortarGridGeometry& mortar, Visitor&& v) const
    {
        std::ranges::for_each(
            mortarToSubDomains_.at(mortarIndexOf_(mortar)),
            [&] (auto sdIdx) { std::visit(v, subDomains_.at(sdIdx)); }
        );
    }

 private:
    std::size_t mortarIndexOf_(const MortarGridGeometry& mortar) const
    {
        auto it = std::ranges::find_if(mortars_, [&] (const auto ptr) { return ptr.get() == &mortar; });
        if (it == mortars_.end())
            DUNE_THROW(Dune::InvalidStateException, "Given mortar grid geometry is not in this decomposition");
        return std::ranges::distance(std::ranges::begin(mortars_), it);
    }

    template<typename GridGeometry>
    std::size_t subDomainIndexOf_(const GridGeometry& subDomain) const
    {
        auto it = std::ranges::find_if(subDomains_, [&] (const auto& sd) {
            return std::visit([&] (const auto& ptr) { return ptr.get() == &subDomain; }, sd);
        });
        if (it == subDomains_.end())
            DUNE_THROW(Dune::InvalidStateException, "Given grid geometry is not in this decomposition");
        return std::ranges::distance(std::ranges::begin(subDomains_), it);
    }

    friend DecompositionFactory<MortarGridGeometry, GridGeometries...>;
    Decomposition() = default;
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
    std::vector<std::vector<std::size_t>> mortarToSubDomains_; // TODO: could use ReservedVector<size_t, 2>?
    std::vector<std::vector<std::size_t>> subDomainsToMortar_;
};

/*!
 * \ingroup MultiDomain
 * \brief Factory for decomposition.
 * \tparam MortarGridGeometry Grid geometry type of the mortar domain.
 * \tparam GridGeometries Types of all subdomain grid geometries that may occur in the decomposition.
 */
template<typename MortarGridGeometry, typename... GridGeometries>
class DecompositionFactory
{
 public:
    //! Insert a mortar domain and return this factory
    DecompositionFactory& withMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { mortars_.push_back(gg); return *this; }

    //! Insert a subdomain and return this factory
    template<typename GridGeometry> requires(std::disjunction_v<std::is_same<std::remove_cvref_t<GridGeometry>, GridGeometries>...>)
    DecompositionFactory& withSubDomain(std::shared_ptr<GridGeometry> gg)
    { subDomains_.push_back(gg); return *this; }

    //! Create a decomposition from all inserted mortars & subdomains
    Decomposition<MortarGridGeometry, GridGeometries...> make() const
    {
        Decomposition<MortarGridGeometry, GridGeometries...> result;
        result.mortars_ = mortars_;
        result.subDomains_ = subDomains_;
        result.mortarToSubDomains_.resize(result.mortars_.size());
        result.subDomainsToMortar_.resize(result.subDomains_.size());
        std::size_t mortarId = 0;
        for (const auto& mortarPtr : result.mortars_)
        {
            std::size_t sdId = 0;
            for (const auto& sd : result.subDomains_)
            {
                const bool intersect = std::visit([&] (const auto& sdPtr) {
                    return makeGlue(*mortarPtr, *sdPtr).size() > 0;
                }, sd);
                if (intersect)
                {
                    result.mortarToSubDomains_[mortarId].push_back(sdId);
                    result.subDomainsToMortar_[sdId].push_back(mortarId);
                }
                sdId++;
            }
            mortarId++;
        }

        return result;
    }

 private:
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
};

} // end namespace Dumux::Mortar

#endif
