// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief A helper function for class member function introspection
 * \note Follows the description by Jean Guegant on
 *       https://jguegant.github.io/blogs/tech/sfinae-introduction.html
 */
#ifndef DUMUX_TYPETRAITS_ISVALID_HH
#define DUMUX_TYPETRAITS_ISVALID_HH

#include <type_traits>

namespace Dumux {

namespace Detail {

// the functor testing an expression for validity
// Expression: can be for example a lambda expression
template <typename Expression>
struct ValidityTestFunctor
{
private:
    // std::declval creates an object of expression
    // the expression, i.e. a lambda expression gets an object as a parameter and does checks on it.
    // so we create also an object of the parameter usiung std::declval
    // if decltype can evaluate the type, i.e. the object parameter is a valid argument for the expression
    // we return std::true_type
    // note: the int is used to give the first overload always precedence
    // note: the last argument in decltype determines the deduced type but all types need to be valid
    template <typename Argument>
    constexpr auto testArgument_(int /* unused int to make this overload the priority choice */) const
    -> decltype(std::declval<Expression>()(std::declval<Argument>()), std::true_type())
    { return std::true_type(); }

    // otherwise we return std::false_type, i.e. this is the fallback
    template <typename Argument>
    constexpr std::false_type testArgument_(...) const
    { return std::false_type(); }

public:
    // the operator () takes the argument we want to use as argument to the test expression
    // note we use the int to prefer the "valid"-result overload of test_ if possible
    template <typename Argument>
    constexpr auto operator() (const Argument& arg) const
    { return testArgument_<Argument>(int()); }

    // check function takes the template argument explicitly
    template <typename Argument>
    constexpr auto check () const
    { return testArgument_<Argument>(int()); }
};

} // end namespace Detail


/*!
 * \ingroup Typetraits
 * \brief A function that creates a test functor to do class member introspection at compile time
 * \return a functor that returns true if the expression is valid with a given type / object
 * Usage:
 * If you want to test if a class has the member function resize(std::size_t) create a test functor
 * \code
 * auto hasResize = isValid([](auto&& c) -> decltype(c.resize(std::size_t(1))) {}; });
 * \endcode
 * \note hasResize can be constexpr in C++17 which allows lambdas in constexpr functions
 * The you can use the test in compile time expressions
 * \code
 * template<class T>
 * auto otherFunc(const T& t)
 * -> typename std::enable_if_t<!decltype(hasResize(myvector))::value, double>
 * { return 4.0; }
 * \endcode
 */
template <typename Expression>
constexpr auto isValid(const Expression& t)
{ return Detail::ValidityTestFunctor<Expression>(); }

} // end namespace Dumux

#endif
