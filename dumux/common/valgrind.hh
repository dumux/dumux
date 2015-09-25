// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Some templates to wrap the valgrind macros
 */
#ifndef DUMUX_VALGRIND_HH
#define DUMUX_VALGRIND_HH

#ifndef HAVE_VALGRIND
// make sure that the HAVE_VALGRIND macro is always defined
#define HAVE_VALGRIND 0
#endif

#if __GNUC__ < 4 || (__GNUC__ == 4  && __GNUC_MINOR__ < 5)
// do not do static_asserts for gcc < 4.5 (semantics changed, so old
// GCCs will complain when using static_assert)
#define static_assert(a, b)

// do not do valgrind client requests for gcc < 4.5 (old GCCs do not
// support anonymous template arguments which results in errors inside
// the BoundaryTypes class)
#define SetUndefined(a) foo();
#define SetDefined(a) foo();
#define CheckDefined(a) foo();
#define SetNoAccess(a) foo();
namespace Valgrind
{
inline void foo() {}; // dummy function
}

#else

#if HAVE_VALGRIND
#include <valgrind/memcheck.h>
#endif // HAVE_VALGRIND

namespace Valgrind
{
/*!
 * \brief Make valgrind complain if the object occupied by an object
 *        is undefined.
 *
 * Please note that this does not check whether the destinations of
 * the object's pointers or references are defined.
 *
 * \tparam T The type of the object which ought to be checked
 *
 * \param value the object which valgrind should check
 */
template <class T>
inline void CheckDefined(const T &value)
{
#if HAVE_VALGRIND
    VALGRIND_CHECK_MEM_IS_DEFINED(&value, sizeof(T));
#endif
}

/*!
 * \brief Make the memory on which an object resides undefined.
 *
 * \tparam T The type of the object which ought to be set to undefined
 *
 * \param value The object which's memory valgrind should be told is undefined
 */
template <class T>
inline void SetUndefined(const T &value)
{
#if HAVE_VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(&value, sizeof(T));
#endif
}

/*!
 * \brief Make the memory on which an object resides defined.
 *
 * \tparam T The type of the object which valgrind should consider as defined
 *
 * \param value The object which's memory valgrind should consider as defined
 */
template <class T>
inline void SetDefined(const T &value)
{
#if HAVE_VALGRIND
    VALGRIND_MAKE_MEM_DEFINED(&value, sizeof(T));
#endif
}

/*!
 * \brief Make valgrind complain if an object's memory is accessed.
 *
 * \tparam T The type of the object which valgrind should complain if accessed
 *
 * \param value The object which's memory valgrind should complain if accessed
 */
template <class T>
inline void SetNoAccess(const T &value)
{
#if HAVE_VALGRIND
    VALGRIND_MAKE_MEM_NOACCESS(&value, sizeof(T));
#endif
}

}

#endif

#endif
