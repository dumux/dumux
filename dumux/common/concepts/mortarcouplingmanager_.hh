// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concept for the coupling manager of a subdomain in a mortar-coupled model
 */
#ifndef DUMUX_CONCEPTS_MORTAR_COUPLING_MANAGER__HH
#define DUMUX_CONCEPTS_MORTAR_COUPLING_MANAGER__HH

#include <concepts>
#include <cstddef>
#include <memory>

#include <dumux/multidomain/mortar/couplingmode.hh>

namespace Dumux::Concept {

/*!
 * \ingroup Concepts
 * \brief The interface a mortar subdomain problem composes: whether a sub-control volume
 *        or a position lies on a mortar trace, the mortar data imposed there and the layout
 *        it is given in, the mode in which it enters the problem, and the state of the
 *        coupled solve. The subdomain solver drives the mutating side.
 */
template<class M>
concept MortarSubDomainCouplingManager = requires(M& manager,
                                                  const M& constManager,
                                                  const typename M::Element& element,
                                                  const typename M::SubControlVolume& scv,
                                                  const typename M::GlobalPosition& globalPos,
                                                  std::size_t mortarId,
                                                  typename M::TraceSolutionVector traceData)
{
    typename M::GridGeometry;
    typename M::Trace;
    typename M::TraceSolutionVector;
    { constManager.gridGeometry() } -> std::convertible_to<const typename M::GridGeometry&>;
    { constManager.couplingMode() } -> std::same_as<Mortar::CouplingMode>;
    { constManager.isHomogeneous() } -> std::convertible_to<bool>;
    { constManager.isFloating() } -> std::convertible_to<bool>;
    { constManager.isCoupled(element, scv) } -> std::convertible_to<bool>;
    { constManager.isCoupledAtPos(globalPos) } -> std::convertible_to<bool>;
    { constManager.isPinnedDof(element, scv) } -> std::convertible_to<bool>;
    { constManager.numTraceCells(mortarId) } -> std::convertible_to<std::size_t>;
    { constManager.traceData(mortarId) } -> std::convertible_to<const typename M::TraceSolutionVector&>;
    { constManager.dataVersion() } -> std::convertible_to<std::size_t>;
    { constManager.traceDataOrder() } -> std::convertible_to<std::size_t>;
    { constManager.numTraceDofs(mortarId) } -> std::convertible_to<std::size_t>;
    constManager.traceAt(element, scv);
    manager.registerTrace(std::shared_ptr<const typename M::Trace>{}, mortarId);
    manager.setTraceDataOrder(std::size_t{});
    manager.setTraceVariables(mortarId, traceData);
    manager.setCouplingMode(Mortar::CouplingMode::essential);
    manager.setHomogeneous(true);
};

} // end namespace Dumux::Concept

#endif
