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
 * \ingroup Common
 * \brief Get only parts of a container or tuple
 */
#ifndef DUMUX_COMMON_PARTIAL_HH
#define DUMUX_COMMON_PARTIAL_HH

#include <tuple>
#include <type_traits>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/common/version.hh>

namespace Dumux {

#if DUNE_VERSION_GT_REV(DUNE_ISTL,2,6,0)
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

#else

// for backwards-compatibility of partial with v2.6.0
template<typename ... Args>
class MultiTypeBlockVectorProxy : public std::tuple<Args...>
{
    typedef std::tuple<Args...> ParentType;
public:
    using ParentType::ParentType;

    typedef MultiTypeBlockVectorProxy<Args...> type;
    typedef double field_type;

    /** \brief Random-access operator
     */
    template< std::size_t index >
    typename std::tuple_element<index, ParentType>::type&
    operator[] ( const std::integral_constant< std::size_t, index > indexVariable )
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }

    /** \brief Const random-access operator
     */
    template< std::size_t index >
    const typename std::tuple_element<index, ParentType>::type&
    operator[] ( const std::integral_constant< std::size_t, index > indexVariable ) const
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }
};

template<class T> struct isMultiTypeBlockVector;

template<class... T>
struct isMultiTypeBlockVector<MultiTypeBlockVectorProxy<T...> > : public std::true_type {};

/*!
 * \brief a function to get a MultiTypeBlockVector with references to some entries of another MultiTypeBlockVector
 * \param v a MultiTypeBlockVector
 * \param indices the indices of the entries that should be referenced
 */
template<class ...Args, std::size_t ...i>
auto partial(Dune::MultiTypeBlockVector<Args...>& v, Dune::index_constant<i>... indices)
{
    return MultiTypeBlockVectorProxy<std::add_lvalue_reference_t<std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(v[indices]...);
}

#endif

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
