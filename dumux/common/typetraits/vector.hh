// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
