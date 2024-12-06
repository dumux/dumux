// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
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

#include <dumux/io/grid/facetgridmanager.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux::Mortar {

template<typename MortarGridGeometry, typename... GridGeometries>
class DecompositionFactory;

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Class to represent a decomposition of a domain into multiple subdomains with mortars on their interfaces.
 * \tparam MortarGridGeometry Grid geometry type of the mortar domain.
 * \tparam GridGeometries Types of all subdomain grid geometries that may occur in the decomposition.
 * \note Construct an instance of this class using a `DecompositionFactory`.
 */
template<typename MortarGridGeometry, typename... GridGeometries>
class Decomposition
{
    using MortarGrid = typename MortarGridGeometry::GridView::Grid;
    using SubDomainTrace = std::variant<std::shared_ptr<FacetGridManager<typename GridGeometries::GridView::Grid, MortarGrid>>...>;

 public:
    std::size_t numberOfMortars() const { return mortars_.size(); }
    std::size_t numberOfSubDomains() const { return subDomains_.size(); }

    //! Return true if the given mortar is part of this decomposition
    bool containsMortar(const MortarGridGeometry& mortar) const {
        return std::ranges::any_of(mortars_, [&] (const auto& ptr) { return ptr.get() == &mortar; });
    }

    //! Return true if the given subdomain is part of this decomposition
    template<typename GridGeometry>
        requires(std::disjunction_v<std::is_same<GridGeometry, GridGeometries>...>)
    bool containsSubDomain(const GridGeometry& gridGeometry) const {
        return std::ranges::any_of(subDomains_, [&] (const auto& sd) {
            return std::visit([&] (const auto& ptr) { return ptr.get() == &gridGeometry; }, sd);
        });
    }

    //! Return the id of the given mortar within this decomposition
    std::size_t id(const MortarGridGeometry& mortar) const
    {
        auto it = std::ranges::find_if(mortars_, [&] (const auto ptr) { return ptr.get() == &mortar; });
        if (it == mortars_.end())
            DUNE_THROW(Dune::InvalidStateException, "Given mortar grid geometry is not in this decomposition");
        return std::ranges::distance(std::ranges::begin(mortars_), it);
    }

    //! Return the id of the given subdomain within this decomposition
    template<typename GridGeometry>
        requires(std::disjunction_v<std::is_same<GridGeometry, GridGeometries>...>)
    std::size_t id(const GridGeometry& subDomain) const
    {
        auto it = std::ranges::find_if(subDomains_, [&] (const auto& sd) {
            return std::visit([&] <typename G> (const std::shared_ptr<const G>& ptr) {
                if constexpr (std::is_same_v<G, GridGeometry>)
                    return ptr.get() == &subDomain;
                return false;
            }, sd);
        });
        if (it == subDomains_.end())
            DUNE_THROW(Dune::InvalidStateException, "Given grid geometry is not in this decomposition");
        return std::ranges::distance(std::ranges::begin(subDomains_), it);
    }

    //! Visit all mortars in this decomposition
    template<std::invocable<const std::shared_ptr<const MortarGridGeometry>&> Visitor>
    void visitMortars(Visitor&& v) const
    { std::ranges::for_each(mortars_, v); }

    //! Visit all subdomains in this decomposition
    template<typename Visitor>
        requires(std::conjunction_v<std::is_invocable<Visitor, const std::shared_ptr<const GridGeometries>&>...>)
    void visitSubDomains(Visitor&& v) const
    { std::ranges::for_each(subDomains_, [&] (const auto& sd) { std::visit(v, sd); }); }

    //! Visit the mortar domains that are coupled to the given subdomain
    template<typename GridGeometry, std::invocable<const std::shared_ptr<const MortarGridGeometry>&> Visitor>
        requires(std::disjunction_v<std::is_same<GridGeometry, GridGeometries>...>)
    void visitCoupledMortarsOf(const GridGeometry& subDomain, Visitor&& v) const
    {
        std::ranges::for_each(
            subDomainsToMortar_.at(id(subDomain)),
            [&] (auto mortarIdx) { v(mortars_.at(mortarIdx)); }
        );
    }

    //! Visit the subdomains that are coupled to the given mortar domain
    template<typename Visitor>
        requires(std::conjunction_v<std::is_invocable<Visitor, const std::shared_ptr<const GridGeometries>&>...>)
    void visitCoupledSubDomainsOf(const MortarGridGeometry& mortar, Visitor&& v) const
    {
        std::ranges::for_each(
            mortarToSubDomains_.at(id(mortar)),
            [&] (auto sdIdx) { std::visit(v, subDomains_.at(sdIdx)); }
        );
    }

    //! Visit the part of the subdomain trace that overlaps with the given mortar
    template<typename GridGeometry, typename Visitor>
        requires(std::disjunction_v<std::is_same<GridGeometry, GridGeometries>...>)
    void visitSubDomainTraceWith(const MortarGridGeometry& mortar, const GridGeometry& subDomain, Visitor&& v) const
    {
        const auto mortarIndex = id(mortar);
        const auto subDomainIndex = id(subDomain);
        const auto& map = subDomainsToMortar_.at(subDomainIndex);
        const auto it = std::ranges::find(map, mortarIndex);
        if (it == std::ranges::end(map))
            DUNE_THROW(Dune::InvalidStateException, "Could not find trace for the given pair of subdomain/mortar.");
        const auto i = std::ranges::distance(std::ranges::begin(map), it);
        std::visit(v, subDomainToMortarTraces_.at(subDomainIndex).at(i));
    }

 private:
    friend DecompositionFactory<MortarGridGeometry, GridGeometries...>;
    Decomposition() = default;
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
    std::vector<std::vector<std::size_t>> mortarToSubDomains_; // TODO: could use ReservedVector<size_t, 2>?
    std::vector<std::vector<std::size_t>> subDomainsToMortar_;
    std::vector<std::vector<SubDomainTrace>> subDomainToMortarTraces_;
};

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Factory for decomposition.
 * \tparam MortarGridGeometry Grid geometry type of the mortar domain.
 * \tparam GridGeometries Types of all subdomain grid geometries that may occur in the decomposition.
 */
template<typename MortarGridGeometry, typename... GridGeometries>
class DecompositionFactory
{
    using MortarGrid = typename MortarGridGeometry::GridView::Grid;

 public:
    //! Insert a mortar domain
    void insertMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { mortars_.push_back(gg); }

    //! Insert a mortar domain and return this factory
    DecompositionFactory& withMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { insertMortar(gg); return *this; }

    //! Insert a subdomain
    template<typename GridGeometry> requires(std::disjunction_v<std::is_same<std::remove_cvref_t<GridGeometry>, GridGeometries>...>)
    void insertSubDomain(std::shared_ptr<GridGeometry> gg)
    { subDomains_.push_back(gg); }

    //! Insert a subdomain and return this factory
    template<typename GridGeometry> requires(std::disjunction_v<std::is_same<std::remove_cvref_t<GridGeometry>, GridGeometries>...>)
    DecompositionFactory& withSubDomain(std::shared_ptr<GridGeometry> gg)
    { insertSubDomain(gg); return *this; }

    //! Create a decomposition from all inserted mortars & subdomains
    Decomposition<MortarGridGeometry, GridGeometries...> make() const
    {
        Decomposition<MortarGridGeometry, GridGeometries...> result;
        result.mortars_ = mortars_;
        result.subDomains_ = subDomains_;
        result.mortarToSubDomains_.resize(result.mortars_.size());
        result.subDomainsToMortar_.resize(result.subDomains_.size());
        result.subDomainToMortarTraces_.resize(result.subDomains_.size());
        std::size_t mortarId = 0;
        for (const auto& mortarPtr : result.mortars_)
        {
            std::size_t sdId = 0;
            unsigned int numNeighbors = 0;
            for (const auto& sd : result.subDomains_)
            {
                const bool intersect = std::visit([&] <typename SD> (const std::shared_ptr<const SD>& sdPtr) {
                    FacetGridManager<typename SD::GridView::Grid, MortarGrid> facetGridManager;
                    facetGridManager.init(sdPtr->gridView().grid(), [&] (const auto&, const auto& is) {
                        return !intersectingEntities(is.geometry(), mortarPtr->boundingBoxTree()).empty();
                    });
                    const bool intersects = facetGridManager.grid().leafGridView().size(0) > 0;
                    if (intersects)
                        result.subDomainToMortarTraces_[sdId].emplace_back(
                            std::make_shared<std::remove_cvref_t<decltype(facetGridManager)>>(std::move(facetGridManager))
                        );
                    return intersects;
                }, sd);
                if (intersect)
                {
                    result.mortarToSubDomains_[mortarId].push_back(sdId);
                    result.subDomainsToMortar_[sdId].push_back(mortarId);
                    numNeighbors++;
                }
                sdId++;
            }
            mortarId++;
            if (numNeighbors != 2)
                DUNE_THROW(Dune::InvalidStateException, "Each mortar is expected to couple to two subdomains");
        }

        return result;
    }

 private:
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
};

} // end namespace Dumux::Mortar

#endif
