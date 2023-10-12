// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Get only parts of a container or tuple
 */
#ifndef DUMUX_COMMON_PARTIAL_HH
#define DUMUX_COMMON_PARTIAL_HH

#include <tuple>
#include <type_traits>

#include <dune/istl/multitypeblockvector.hh>

namespace Dumux {

/*!
 * \brief a function to get a MultiTypeBlockVector with references to some entries of another MultiTypeBlockVector
 * \param v a MultiTypeBlockVector
 * \param indices the indices of the entries that should be referenced
 */
template<class ...Args, std::size_t ...i>
auto partial(Dune::MultiTypeBlockVector<Args...>& v, Dune::index_constant<i>... indices)
{
    return Dune::MultiTypeBlockVector<std::add_lvalue_reference_t<std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(v[indices]...);
}

/*!
 * \brief a function to get a MultiTypeBlockVector with const references to some entries of another MultiTypeBlockVector
 * \param v a MultiTypeBlockVector
 * \param indices the indices of the entries that should be referenced
 */
template<class ...Args, std::size_t ...i>
auto partial(const Dune::MultiTypeBlockVector<Args...>& v, Dune::index_constant<i>... indices)
{
    return Dune::MultiTypeBlockVector<std::add_lvalue_reference_t<const std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(v[indices]...);
}

/*!
 * \brief a function to get a tuple with references to some entries of another tuple
 * \param v a tuple
 * \param indices a tuple of indices of the entries that should be referenced
 */
template<class ...Args, std::size_t ...i>
auto partial(std::tuple<Args...>& v, Dune::index_constant<i>... indices)
{
    return std::tuple<std::add_lvalue_reference_t<std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(std::get<indices>(v)...);
}

/*!
 * \brief a function to get a tuple with const references to some entries of another tuple
 * \param v a tuple
 * \param indices a tuple of indices of the entries that should be referenced
 */
template<class ...Args, std::size_t ...i>
auto partial(const std::tuple<Args...>& v, Dune::index_constant<i>... indices)
{
    return std::tuple<std::add_lvalue_reference_t<const std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(std::get<indices>(v)...);
}

/*!
 * \brief a function to get a MultiTypeBlockVector with references to some entries of another MultiTypeBlockVector
 * \param t an std::tuple or Dune::MultiTypeBlockVector
 * \param indices a tuple of indices of the entries that should be referenced
 */
template<class T, std::size_t ...i>
auto partial(T& t, std::tuple<Dune::index_constant<i>...> indices)
{
    return partial(t, Dune::index_constant<i>{}...);
}

} // end namespace Dumux

#endif
