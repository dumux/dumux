// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits to be used for detecting new interfaces related to local dofs
 */
#ifndef DUMUX_TYPETRAITS_LOCALDOFS__HH
#define DUMUX_TYPETRAITS_LOCALDOFS__HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>
#include <dune/common/rangeutilities.hh>

namespace Dumux::Detail::LocalDofs {

//! helper struct detecting if a fvElementGeometry object has a numLocalDofs() function
template<class Imp>
using NumLocalDofsDetector = decltype(
    std::declval<Imp>().numLocalDofs()
);

template<typename ElementDiscretization>
constexpr int numLocalDofs(const ElementDiscretization& fvGeometry)
{
    if constexpr (Dune::Std::is_detected<NumLocalDofsDetector, ElementDiscretization>::value)
        return fvGeometry.numLocalDofs();
    else
        return fvGeometry.numScv();
}

//! helper struct detecting if a fvElementGeometry object has a nonCVLocalDofs() function
template<class Imp>
using NonCVLocalDofsDetector = decltype(
    nonCVLocalDofs(std::declval<Imp>())
);

template<class Imp>
constexpr inline bool hasNonCVLocalDofsInterface()
{ return Dune::Std::is_detected<NonCVLocalDofsDetector, Imp>::value; }

//! helper struct detecting if a fvElementGeometry object has maxNumElementDofs
template<class Imp>
using MaxNumElementDofs = decltype( Imp::maxNumElementDofs );

template<typename ElementDiscretization>
constexpr int maxNumLocalDofs()
{
    if constexpr (Dune::Std::is_detected<MaxNumElementDofs, ElementDiscretization>::value)
        return ElementDiscretization::maxNumElementDofs;
    else
        return ElementDiscretization::maxNumElementScvs;
}

//! helper struct detecting if a class has a localDofIndex() function
template<class Imp>
using LocalDofIndexDetector = decltype(
    std::declval<Imp>().localDofIndex()
);

template<class ScvOrLocalDof>
inline auto index(const ScvOrLocalDof& scvOrLocalDof)
{
    if constexpr (Dune::Std::is_detected<LocalDofIndexDetector, ScvOrLocalDof>::value)
        return scvOrLocalDof.localDofIndex();
    else
        return scvOrLocalDof.index();
}

} // end namespace Dumux::Detail::LocalDofs

#endif
