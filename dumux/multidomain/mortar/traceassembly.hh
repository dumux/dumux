// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Assembly of the trace a subdomain hands back to the mortar domain.
 *
 * A trace is the mean of a quantity over each trace cell. It is assembled from the
 * sub-entities of the subdomain that lie on the cell, each contributing the integral of the
 * quantity over itself, and the sum is divided by the cell's measure. An integrand names the
 * kind of sub-entity it is evaluated on and returns that integral; the integrands of the
 * physical models live next to those models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_TRACE_ASSEMBLY_HH
#define DUMUX_MULTIDOMAIN_MORTAR_TRACE_ASSEMBLY_HH

#include <cstddef>
#include <type_traits>
#include <utility>

#include <dumux/common/concepts/mortarcouplingmanager_.hh>
#include <dumux/flux/tracefields.hh>

#include "couplingmanager.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief The integral of a quantity over a sub-entity of the given kind, evaluated on the
 *        subdomain's element context bound to the sub-entity's element.
 */
template<TraceEntity e, typename F>
struct TraceIntegrand
{
    static constexpr TraceEntity entity = e;

    //! The integral over the given sub-entity, in the context bound to its element
    template<typename Context, typename SubEntity>
    auto operator()(const Context& context, const SubEntity& subEntity) const
    { return integrand(context, subEntity); }

    F integrand;
};

//! Tag a callable `(context, subEntity) -> integral` with the kind of sub-entity it integrates over
template<TraceEntity entity, typename F>
TraceIntegrand<entity, std::decay_t<F>> traceIntegrand(F&& integrand)
{ return {std::forward<F>(integrand)}; }

/*!
 * \ingroup MortarCoupling
 * \brief The trace of a quantity over the cells of the trace shared with the given mortar
 *        domain, as the mortar projectors consume it.
 *
 * Which quantity to integrate is the modelling decision of the conjugate trace: the
 * quantity the mortar residual tests against the mortar basis, conjugate to the data
 * imposed on the trace. Evaluating it with the quadrature the subdomain's own boundary
 * conditions use makes the trace read back the adjoint of the imposition.
 */
template<Concept::MortarSubDomainCouplingManager CouplingManager, typename Problem,
         typename GridVariables, typename SolutionVector, TraceEntity entity, typename F>
auto assembleTrace(const CouplingManager& couplingManager,
                   std::size_t mortarId,
                   const Problem& problem,
                   const GridVariables& gridVariables,
                   const SolutionVector& x,
                   TraceIntegrand<entity, F> integrand)
{
    return traceValues<entity>(
        couplingManager, mortarId,
        makeNonOwningFluxField(problem, gridVariables, x, std::move(integrand))
    );
}

} // end namespace Dumux::Mortar

#endif
