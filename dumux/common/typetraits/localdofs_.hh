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

#include <dune/common/std/type_traits.hh>
#include <dune/common/rangeutilities.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Detail::LocalDofs {

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

//! helper struct detecting if a fvElementGeometry object defines its own local dof type
template<class FVG>
using SpecifiesLocalDof = typename FVG::LocalDof;

template<class FVG>
using LocalDof_t = Dune::Std::detected_or_t<
    Dumux::CVFE::LocalDof<typename IndexTraits<typename FVG::GridGeometry::GridView>::LocalIndex,
                          typename IndexTraits<typename FVG::GridGeometry::GridView>::GridIndex >,
    SpecifiesLocalDof,
    FVG
>;

} // end namespace Dumux::Detail::LocalDofs

#endif
