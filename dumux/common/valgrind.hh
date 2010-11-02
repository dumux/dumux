// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
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
