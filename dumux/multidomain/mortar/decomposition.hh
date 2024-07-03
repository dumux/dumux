// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Class that represents a domain decomposition.
 *        Contains the subdomain & mortar grid geometries and connectivity information between them.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH
#define DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH

#include <tuple>
#include <array>
#include <vector>
#include <variant>
#include <utility>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/multidomain/glue.hh>
#include <dumux/multidomain/fvgridgeometry.hh>

namespace Dumux::Mortar {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief Class that represents the decomposition of a domain.
 * \tparam SubDomainGridGeometries A `MultiDomainFVGridGeometry` containing the subdomain grid geometries.
 ' \tparam MortarGridGeometries A `MultiDomainFVGridGeometry` containing the mortar grid geometries.
 */
template<typename SubDomainGridGeometries, typename MortarGridGeometries>
class Decomposition
{
    template<typename T>
    struct DomainIndexVariant;

    template<std::size_t... i>
    struct DomainIndexVariant<std::index_sequence<i...>>
    : std::type_identity<std::variant<Dune::index_constant<i>...>>
    {};

    // class to store either a reference or an instance of std::remove_reference_t<T>
    template<typename T>
    class Storage
    {
        static constexpr bool isConst = std::is_const_v<std::remove_reference_t<T>>;
        static constexpr bool isReference = std::is_lvalue_reference_v<T>;

    public:
        using Type = std::remove_cvref_t<T>;

        Storage(T& t) requires(isReference) : stored_{t} {}
        Storage(Type&& t) requires(!isReference) : stored_{std::move(t)} {}

        const Type& get() const { return stored_; }
        Type& get() requires(!isConst) { return stored_; }

    private:
        std::conditional_t<isReference, T, std::remove_cvref_t<T>> stored_;
    };

    using SubDomainStorage = Storage<SubDomainGridGeometries>;
    using MortarStorage = Storage<MortarGridGeometries>;

public:
    static constexpr std::size_t numSubDomains = SubDomainStorage::Type::size;
    static constexpr std::size_t numMortars = MortarStorage::Type::size;

    template<std::size_t i>
    using SubDomainGridGeometry = typename SubDomainStorage::Type::template Type<i>;

    template<std::size_t i>
    using MortarGridGeometry = typename MortarStorage::Type::template Type<i>;

    template<typename SDGG, typename MGG>
        requires(std::is_same_v<std::remove_cvref_t<SDGG>, typename SubDomainStorage::Type> and
                 std::is_same_v<std::remove_cvref_t<MGG>, typename MortarStorage::Type>)
    explicit Decomposition(SDGG&& subdomainGGs, MGG&& mortarGGs)
    : subdomainGridGeometries_{std::forward<SDGG>(subdomainGGs)}
    , mortarGridGeometries_{std::forward<MGG>(mortarGGs)} {
        Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<numSubDomains>{}), [&] (auto&& sdId) {
            Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<numMortars>{}), [&] (auto&& mId) {
                const auto glue = makeGlue(subDomainGridGeometry(sdId), mortarGridGeometry(mId));

                if (glue.size() > 0)
                {
                    std::cout << "Found intersections between subdomain " << sdId << " and mortar " << mId << std::endl;
                    subDomainToMortarIds_[sdId].push_back({mId});
                    mortarToSubDomainIds_[mId].push_back({sdId});
                    if (mortarToSubDomainIds_[mId].size() > 2)
                        DUNE_THROW(Dune::InvalidStateException, "Each mortar domain must only be between two subdomains");
                }
            });
        });
    }

    template<std::size_t i>
    const SubDomainGridGeometry<i>& subDomainGridGeometry(const Dune::index_constant<i>& id) const
    { return subdomainGridGeometries_.get()[id]; }

    template<std::size_t i>
    SubDomainGridGeometry<i>& subDomainGridGeometry(const Dune::index_constant<i>& id)
    { return subdomainGridGeometries_.get()[id]; }

    template<std::size_t i>
    const MortarGridGeometry<i>& mortarGridGeometry(const Dune::index_constant<i>& id) const
    { return mortarGridGeometries_.get()[id]; }

    template<std::size_t i>
    MortarGridGeometry<i>& mortarGridGeometry(const Dune::index_constant<i>& id)
    { return mortarGridGeometries_.get()[id]; }

private:
    using SubDomainIndexVariant = typename DomainIndexVariant<std::make_index_sequence<numSubDomains>>::type;
    using MortarIndexVariant = typename DomainIndexVariant<std::make_index_sequence<numMortars>>::type;

    Storage<SubDomainGridGeometries> subdomainGridGeometries_;
    Storage<MortarGridGeometries> mortarGridGeometries_;

    std::array<std::vector<SubDomainIndexVariant>, numSubDomains> subDomainToMortarIds_;
    std::array<std::vector<MortarIndexVariant>, numMortars> mortarToSubDomainIds_;
};

template<typename SDGG, typename MGG>
Decomposition(SDGG&&, MGG&&) -> Decomposition<SDGG, MGG>;

} // end namespace Dumux::Mortar

#endif
