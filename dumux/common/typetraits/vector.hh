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
 * \ingroup Typetraits
 * \brief Type traits to be used with vector types
 */
#ifndef DUMUX_TYPETRAITS_VECTOR_HH
#define DUMUX_TYPETRAITS_VECTOR_HH

#include <type_traits>
#include <dune/istl/multitypeblockvector.hh>


namespace Dumux {

    //! Helper type to determine whether a given type is a Dune::MultiTypeBlockVector
    template<class T> struct isMultiTypeBlockVector : public std::false_type {};

    //! Helper type to determine whether a given type is a Dune::MultiTypeBlockVector
    template<class... T>
    struct isMultiTypeBlockVector<Dune::MultiTypeBlockVector<T...> > : public std::true_type {};

} // end namespace Dumux

#endif
