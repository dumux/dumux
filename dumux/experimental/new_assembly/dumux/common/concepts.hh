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
 * \brief Basic concepts.
 */
#ifndef DUMUX_COMMON_CONCEPTS_HH
#define DUMUX_COMMON_CONCEPTS_HH

#include <concepts>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/traits.hh>
#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>

namespace Dumux::Concepts {

//! Concept for arithmetic scalar types
template<typename T>
concept Arithmetic
    = std::integral<T>
    or std::floating_point<T>
    or Traits::IsArithmetic<T>::value;


//! Concept for (two) addable types
template<typename T, typename V = T>
concept Addable = requires(const T& t, const V& v) {
    { t + v };
};


//! Concept for in-place addable types
template<typename T, typename V = T>
concept AddAssignable = requires(T& t, const V& v) {
    { t += v };
};


//! Concept for in-place scalable types
template<typename T, typename Scalar = double>
concept ScaleAssignable = requires(T& t, const Scalar& v) {
    { t *= v };
};


//! Concept for types that allow indexing with `operator[]`
template<typename T, typename I = std::size_t>
concept Indexable = isIndexable<T, I>;


//! Concept for multi-dimensional arrays
template<typename T, int expectedIndexDepth = -1, typename Index = std::size_t>
concept MDArray = Indexable<T, Index> and (
    expectedIndexDepth == -1 or
    indexDepth<T> == expectedIndexDepth
);


//! Concept for types that are tagged as views.
template<typename T>
concept View = Traits::IsView<T>::value;


//! Concept for types that we accept for defining the sizes of containers
template<typename T>
concept Size = std::integral<T> or std::same_as<T, DynamicSize>;


//! Functions export the domain and range types and allow for evaluation at a given point.
template<typename T>
concept Function = requires {
    typename T::Domain;
    typename T::Range;
} and requires(T& t, const typename T::Domain& evalPoint) {
    { t.evaluateAt(evalPoint) } -> std::convertible_to<typename T::Range>;
};


//! Concept of a linearization, exposing the derivative and function value at an evaluation point.
template<typename T>
concept Linearization = requires(const T& t) {
    typename T::Derivative;
    typename T::Value;
    { t.derivative() } -> std::convertible_to<typename T::Derivative>;
    { t.value() } -> std::convertible_to<typename T::Value>;
};


//! Concept for functions that allow for linearization
template<typename T>
concept LinearizableFunction
    = Function<T>
    and requires {
        typename T::Linearization;
        Linearization<typename T::Linearization>;
    } and requires(T& t, const typename T::Domain& evalPoint) {
        { t.linearizeAt(evalPoint) } -> std::convertible_to<typename T::Linearization>;
};


//! Concept for time levels, exposing information on the current time
//! and (potentially) on an ongoing time integration step
template<typename T>
concept TimeLevel = requires(const T& t) {
    typename T::Scalar;

    std::constructible_from<T, typename T::Scalar>; // ctor w/o time step info
    std::constructible_from<T,
                            typename T::Scalar,
                            typename T::Scalar,
                            typename T::Scalar>; // ctor with time step info

    { t.current() } -> std::convertible_to<typename T::Scalar>;
    { t.previous() } -> std::convertible_to<typename T::Scalar>;
    { t.timeStepFraction() } -> std::convertible_to<typename T::Scalar>;
};

} // namespace Dumux::Concepts

#endif
