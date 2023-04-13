// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Multidomain wrapper for multiple grid geometries
 */
#ifndef DUMUX_MULTIDOMAIN_FV_GRIDGEOMETRY_HH
#define DUMUX_MULTIDOMAIN_FV_GRIDGEOMETRY_HH

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
     * \brief Construct grid geometries for all subdomains
     * \param args a list of arguments to pass to the constructors
     *
     * The number of arguments has to match the number of subdomains.
     * In case a constructor needs multiple arguments, they have to be wrapped in a std::tuple.
     * Use std::make_tuple and possible wrap arguments using std::ref / std::cref or use std::forward_as_tuple.
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
    : gridGeometries_(std::move(ggTuple))
    {}

    /*!
     * \brief Update all grid geometries (do this e.g. after grid adaption)
     * \param args a list of arguments to pass to the update functions
     *
     * The number of arguments has to match the number of subdomains.
     * In case the update function needs multiple arguments, they have to be wrapped in a std::tuple.
     * Use std::make_tuple and possible wrap arguments using std::ref / std::cref or use std::forward_as_tuple.
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

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i>) const
    { return *std::get<i>(gridGeometries_); }

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i>)
    { return *std::get<i>(gridGeometries_); }

    ///! access the grid geometry pointer for domain with index i
    template<std::size_t i>
    const PtrType<i>& get(Dune::index_constant<i> id = Dune::index_constant<i>{}) const
    { return std::get<i>(gridGeometries_); }

    ///! access the the grid geometry pointer for domain with index i
    template<std::size_t i>
    PtrType<i>& get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return std::get<i>(gridGeometries_); }

    /*!
     * \brief Access the underlying tuple representation
     */
    TupleType& asTuple()
    { return gridGeometries_; }

    /*!
     * \brief Access the underlying tuple representation
     */
    const TupleType& asTuple() const
    { return gridGeometries_; }

private:

    //! a tuple of points to all grid variables
    TupleType gridGeometries_;
};

} // end namespace Dumux

#endif
