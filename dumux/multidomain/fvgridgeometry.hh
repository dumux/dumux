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
 * \ingroup MultiDomain
 * \brief Multidomain wrapper for multiple grid geometries
 */
#ifndef DUMUX_MULTIDOMAIN_FVGRIDGEOMETRY_HH
#define DUMUX_MULTIDOMAIN_FVGRIDGEOMETRY_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

namespace Dumux {

namespace Multidomain::Detail {

template<class T>
struct IsStdTuple_
: public std::false_type {};

template<class... Args>
struct IsStdTuple_<std::tuple<Args...>>
: public std::true_type {};

template<class T>
inline constexpr bool isStdTuple = IsStdTuple_<T>::value;

} // end namespace Multidomain::Detail

/*!
 * \ingroup MultiDomain
 * \brief A multidomain wrapper for multiple grid geometries
 * \tparam MDTraits The multidomain traits
 */
template<class MDTraits>
class MultiDomainFVGridGeometry
{
    static constexpr std::size_t numSubDomains = MDTraits::numSubDomains;

    // unwrap the tuple and pass its elements as arguments to the constructor of the ith element
    template<std::size_t i, class Tuple, size_t... Is>
    void constructFromTupleOfArgs_(Tuple&& t, std::index_sequence<Is...>)
    {
        std::get<i>(gridGeometries_) = std::make_shared<Type<i>>(std::get<Is>(std::forward<Tuple>(t))...);
    }

    // construct the ith element in this multidomain wrapper
    template<std::size_t i, class Arg>
    void construct_(Arg&& arg)
    {
        using ArgT = std::decay_t<Arg>;
        // use perfect forwarding of the argument(s) in both cases
        if constexpr (Multidomain::Detail::template isStdTuple<ArgT>)
            constructFromTupleOfArgs_<i>(std::forward<Arg>(arg), std::make_index_sequence<std::tuple_size_v<ArgT>>{});
        else
            std::get<i>(gridGeometries_) = std::make_shared<Type<i>>(std::forward<Arg>(arg));
    }

    // unwrap the tuple and pass its elements as arguments to the update function of the ith element
    template<std::size_t i, class Tuple, size_t... Is>
    void updateWithTupleOfArgs_(Tuple&& t, std::index_sequence<Is...>)
    {
        std::get<i>(gridGeometries_)->update(std::get<Is>(std::forward<Tuple>(t))...);
    }

    // update the ith element in this multidomain wrapper
    template<std::size_t i, class Arg>
    void update_(Arg&& arg)
    {
        using ArgT = std::decay_t<Arg>;
        // use perfect forwarding of the argument(s) in both cases
        if constexpr (Multidomain::Detail::template isStdTuple<ArgT>)
            updateWithTupleOfArgs_<i>(std::forward<Arg>(arg), std::make_index_sequence<std::tuple_size_v<ArgT>>{});
        else
            std::get<i>(gridGeometries_)->update(std::forward<Arg>(arg));
    }

public:
    //! export base types of the stored type
    template<std::size_t i>
    using Type = typename MDTraits::template SubDomain<i>::GridGeometry;

    //! export pointer types the stored type
    template<std::size_t i>
    using PtrType = std::shared_ptr<Type<i>>;

    //! export type of tuple of pointers
    using TupleType = typename MDTraits::template Tuple<PtrType>;

    /*!
     * \brief The default constructor
     */
    [[deprecated("Will be removed after release 3.5. Use variadic constructor!")]]
    MultiDomainFVGridGeometry() = default;

    /*!
     * \brief Construct grid geometries for all subdomains
     * \param args a list of arguments to pass to the constructors
     *
     * The number of arguments has to match the number of subdomains.
     * In case a constructor needs multiple arguments, they have to be wrapped in a std::tuple.
     * Use std::make_tuple and possible wrap arguments using std::ref/std::cref or use std::forward_as_tuple.
     * If an argument is a tuple, it will be unpacked and its members will be passed to the constructor.
     * In the corner case where you need to pass a tuple to the constructor,
     * you therefore need to additionally wrap the tuple in a tuple before passing.
     */
    template<typename... Args>
    MultiDomainFVGridGeometry(Args&&... args)
    {
        static_assert(numSubDomains == sizeof...(Args), "Number of arguments has to match number of subdomains");

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [this, t = std::forward_as_tuple(args...)](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            this->construct_<i>(std::get<i>(t));
        });
    }

    /*!
     * \brief Construct wrapper from a tuple of grid geometries
     * \param ggTuple a tuple of shared_ptrs to the grid geometries
     */
    MultiDomainFVGridGeometry(TupleType ggTuple)
    : gridGeometries_(ggTuple)
    {}

    /*!
     * \brief Update all grid geometries (do this e.g. after grid adaption)
     * \param args a list of arguments to pass to the update functions
     *
     * The number of arguments has to match the number of subdomains.
     * In case the update function needs multiple arguments, they have to be wrapped in a std::tuple.
     * Use std::make_tuple and possible wrap arguments using std::ref/std::cref or use std::forward_as_tuple.
     * If an argument is a tuple, it will be unpacked and its members will be passed to the constructor.
     * In the corner case where you need to pass a tuple to the constructor,
     * you therefore need to additionally wrap the tuple in a tuple before passing.
     */
    template<typename... Args>
    void update(Args&&... args)
    {
        static_assert(numSubDomains == sizeof...(Args), "Number of arguments has to match number of subdomains");

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [this, t = std::forward_as_tuple(args...)](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            this->update_<i>(std::get<i>(t));
        });
    }

    /*!
     * \brief Update all grid geometries (do this again after grid adaption)
     */
    [[deprecated("Will be removed after release 3.5. Use variadic update!")]]
    void update()
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridGeometries_, id)->update();
        });
    }

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i>) const
    { return *std::get<i>(gridGeometries_); }

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i>)
    { return *std::get<i>(gridGeometries_); }

    ///! create a copy of the grid geometry pointer for domain with index i
    template<std::size_t i>
    PtrType<i> get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return std::get<i>(gridGeometries_); }

    //! set the pointer for sub domain i
    template<std::size_t i>
    [[deprecated("Will be removed after release 3.5. Use one of the constructors instead.")]]
    void set(PtrType<i> p, Dune::index_constant<i> id = Dune::index_constant<i>{})
    { std::get<i>(gridGeometries_) = p; }

    /*!
     * \brief return the grid variables tuple we are wrapping
     * \note the copy is not expensive since it is a tuple of shared pointers
     */
    TupleType getTuple()
    { return gridGeometries_; }

private:

    //! a tuple of pointes to all grid variables
    TupleType gridGeometries_;
};

} // end namespace Dumux

#endif
