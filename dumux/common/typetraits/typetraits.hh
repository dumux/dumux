// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits.
 */
#ifndef DUMUX_TYPE_TRAITS_HH
#define DUMUX_TYPE_TRAITS_HH

#include <type_traits>

namespace Dumux {

/*!
 * \brief Template which always yields a false value
 * \tparam T Some type.
 */
template<typename T>
struct AlwaysFalse : public std::false_type {};

/*!
 * \brief Function that performs no operation.
 */
inline constexpr auto noop = [] (auto...) {};
using Noop = decltype(noop);

/*!
 * \brief Helper template to select type T if it is not void
 *        or fall back to the given default type otherwise.
 */
template<typename Default, typename T>
using NonVoidOr = std::conditional_t<!std::is_void_v<T>, T, Default>;

/*!
 * \brief Type trait to check if a given type is contained in the provided list of types
 */
template<typename T, typename... Ts>
struct IsAnyOf : std::bool_constant<std::disjunction_v<std::is_same<T, Ts>...>> {};

#ifndef DOXYGEN
namespace Detail {

template<typename... T>
struct TypeList {};

// primary template / entry point
template<typename T, typename... Ts>
struct UniqueTypes {
    using type = std::conditional_t<
        IsAnyOf<T, Ts...>::value,
        typename UniqueTypes<Ts...>::type,
        typename UniqueTypes<TypeList<T>, Ts...>::type
    >;
};

// primary template specialization for empty packs
template<typename T>
struct UniqueTypes<T> : std::type_identity<TypeList<T>> {};

template<typename... Ts, typename T, typename... Rest>
struct UniqueTypes<TypeList<Ts...>, T, Rest...> {
    using type = std::conditional_t<
        IsAnyOf<T, Ts...>::value,
        typename UniqueTypes<TypeList<Ts...>, Rest...>::type,
        typename UniqueTypes<TypeList<Ts..., T>, Rest...>::type
    >;
};

// specizalization for stopping the recursion
template<typename... Ts>
struct UniqueTypes<TypeList<Ts...>> : std::type_identity<TypeList<Ts...>> {};

template<template<typename...> typename T, typename Args>
struct WithTemplateParameters;
template<template<typename...> typename T, typename... Args>
struct WithTemplateParameters<T, TypeList<Args...>> : std::type_identity<T<Args...>> {};

}  // namespace Detail
#endif  // DOXYGEN

/*!
 * \brief Type trait to expose a variadically templated type with unique template parameters.
 */
template<template<typename...> typename T, typename... Args>
struct WithUniqueTemplateParameters : Detail::WithTemplateParameters<T, typename Detail::UniqueTypes<Args...>::type> {};

} // end namespace Dumux

#endif
