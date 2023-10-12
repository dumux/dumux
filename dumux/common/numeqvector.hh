// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief A helper to deduce a vector with the same size as numbers of equations
 */
#ifndef DUMUX_COMMON_NUMEQVECTOR_HH
#define DUMUX_COMMON_NUMEQVECTOR_HH

#include <cstddef>

namespace Dumux {

template<class PrimaryVariables>
struct NumEqVectorTraits
{
    static constexpr std::size_t numEq = PrimaryVariables::size();
    using type = PrimaryVariables;
};

/*!
 * \ingroup Core
 * \brief A vector with the same size as numbers of equations
 * This is the default implementation and has to be specialized for
 * all custom primary variable vector types
 * \note This is based on the primary variables concept
 */
template<class PrimaryVariables>
using NumEqVector = typename NumEqVectorTraits<PrimaryVariables>::type;

} // namespace Dumux

#endif
