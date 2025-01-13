// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits to be used for detecting new interfaces related to local dofs
 */
#ifndef DUMUX_TYPETRAITS_LOCALDOFS_HH
#define DUMUX_TYPETRAITS_LOCALDOFS_HH

#include <dune/common/std/type_traits.hh>
#include <dune/common/rangeutilities.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Detail {

//! helper struct detecting if a fvElementGeometry object has a numLocalDofs() function
template<class Imp>
using NumLocalDofsDetector = decltype(
    std::declval<Imp>().numLocalDofs()
);

template<typename FVElementGeometry>
constexpr int numLocalDofs(const FVElementGeometry& fvGeometry)
{
    if constexpr (Dune::Std::is_detected<NumLocalDofsDetector, FVElementGeometry>::value)
        return fvGeometry.numLocalDofs();
    else
        return fvGeometry.numScv();
}

//! helper struct detecting if a fvElementGeometry object has maxNumElementDofs
template<class Imp>
using MaxNumElementDofs = decltype( Imp::maxNumElementDofs );

template<typename FVElementGeometry>
constexpr int maxNumLocalDofs()
{
    if constexpr (Dune::Std::is_detected<MaxNumElementDofs, FVElementGeometry>::value)
        return FVElementGeometry::maxNumElementDofs;
    else
        return FVElementGeometry::maxNumElementScvs;
}

} // end namespace Dumux::Detail

#endif
