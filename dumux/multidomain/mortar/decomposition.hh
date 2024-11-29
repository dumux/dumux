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
#include <variant>
#include <memory>
#include <vector>

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


 private:
    friend DecompositionFactory<MortarGridGeometry, GridGeometries...>;
    Decomposition() = default;
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
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
    DecompositionFactory& withMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { mortars_.push_back(gg); return *this; }

    template<typename GridGeometry> requires(std::disjunction_v<std::is_same<std::remove_cvref_t<GridGeometry>, GridGeometries>...>)
    DecompositionFactory& withSubDomain(std::shared_ptr<GridGeometry> gg)
    { subDomains_.push_back(gg); return *this; }

    Decomposition<MortarGridGeometry, GridGeometries...> make() const
    {
        Decomposition<MortarGridGeometry, GridGeometries...> result;
        result.mortars_ = mortars_;
        result.subDomains_ = subDomains_;
        // TODO: set up mappings
        return result;
    }

 private:
    std::vector<std::shared_ptr<const MortarGridGeometry>> mortars_;
    std::vector<std::variant<std::shared_ptr<const GridGeometries>...>> subDomains_;
};

} // end namespace Dumux::Mortar

#endif
