// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Helper functions to be used for memory allocation operations on matrices & vectors
 *        that can be called for both static and dynamic types. This is useful wherever it
 *        is not always known if the matrix & vector types at hand are static or dynamic.
 */
#ifndef DUMUX_COMMON_MATRIX_VECTOR_HELPER_HH
#define DUMUX_COMMON_MATRIX_VECTOR_HELPER_HH

#include <type_traits>

namespace Dumux {

//! Determines whether or not a matrix has a resize() function
template<typename M>
struct matrix_has_resize_method
{
private:
    typedef std::true_type yes;
    typedef std::false_type no;

    // resize function is called with two indices for matrices
    template<typename U> static auto test(int) -> decltype(std::declval<U>().resize(0, 0), yes());
    template<typename> static no test(...);

public:
    static constexpr bool value = std::is_same<decltype(test<M>(0)), yes>::value;
};

//! determines whether or not a vector has a resize() function
template<typename V>
struct vector_has_resize_method
{
private:
    typedef std::true_type yes;
    typedef std::false_type no;

    // resize function is called with one index for vectors
    template<typename U> static auto test(int) -> decltype(std::declval<U>().resize(0), yes());
    template<typename> static no test(...);

public:
    static constexpr bool value = std::is_same<decltype(test<V>(0)), yes>::value;
};

//! resizes a matrix to the given sizes (specialization for dynamic matrix type)
template< class Matrix,
          class size_type,
          std::enable_if_t<matrix_has_resize_method<Matrix>::value, int> = 0 >
void resizeMatrix(Matrix& M, size_type rows, size_type cols)
{
    M.resize(rows, cols);
}

//! resizes a matrix to the given sizes (specialization for static matrix type - do nothing)
template< class Matrix,
          class size_type,
          std::enable_if_t<!matrix_has_resize_method<Matrix>::value, int> = 0 >
void resizeMatrix(Matrix& M, size_type rows, size_type cols) {}

//! resizes a vector to the given size (specialization for dynamic matrix type)
template< class Vector,
          class size_type,
          std::enable_if_t<vector_has_resize_method<Vector>::value, int> = 0 >
void resizeVector(Vector& v, size_type size)
{
    v.resize(size);
}

//! resizes a vector to the given size (specialization for static vector type - do nothing)
template< class Vector,
          class size_type,
          std::enable_if_t<!vector_has_resize_method<Vector>::value, int> = 0 >
void resizeVector(Vector& v, size_type rows) {}

} // end namespace Dumux

#endif
