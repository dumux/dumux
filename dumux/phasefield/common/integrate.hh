// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseField
 * \brief Integrate a scalar function of the interpolated field(s) over the
 *        whole domain of a hybrid CVFE (PQ2/PQ3) discretization.
 */
#ifndef DUMUX_PHASEFIELD_COMMON_INTEGRATE_HH
#define DUMUX_PHASEFIELD_COMMON_INTEGRATE_HH

#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux::PhaseField {

/*!
 * \ingroup PhaseField
 * \brief Integrate a scalar function of the interpolated field(s) over the
 *        whole domain, e.g. for mass-conservation or energy-dissipation checks.
 *
 * \param gridGeometry The grid geometry
 * \param gridVariablesCache The (current or previous) grid variables cache
 * \param sol The solution vector matching `gridVariablesCache`
 * \param valueAtIp A callable `(fvGeometry, elemVars, ipCache) -> Scalar`
 *        evaluating the integrand at an interpolation point, typically by
 *        constructing the model's own flux-function-context wrapper.
 */
template<class GridGeometry, class GridVariablesCache, class SolutionVector, class ValueAtIp>
auto integrateField(const GridGeometry& gridGeometry,
                    const GridVariablesCache& gridVariablesCache,
                    const SolutionVector& sol,
                    const ValueAtIp& valueAtIp)
{
    using Scalar = typename SolutionVector::block_type::value_type;

    auto fvGeometry = localView(gridGeometry);
    auto elemVars = localView(gridVariablesCache);

    Scalar integral = 0.0;
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        fvGeometry.bind(element);
        elemVars.bind(element, fvGeometry, sol);

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
        {
            const auto& ipData = qpData.ipData();
            const auto& ipCache = cache(elemVars, ipData);
            integral += qpData.weight()*valueAtIp(fvGeometry, elemVars, ipCache);
        }
    }

    return integral;
}

/*!
 * \ingroup PhaseField
 * \brief Integrate a scalar function of the control-volume (vertex) dofs
 *        *only* over the whole domain, weighted exactly as the model's own
 *        `storageIntegral` weights them (`Extrusion::volume(fvGeometry, scv)`).
 *
 * Unlike integrateField() (which integrates the full, consistent hybrid
 * CVFE field, including non-CV/edge/center dofs), this deliberately
 * restricts itself to control-volume dofs, because the hybrid scheme's flux
 * discretization at those dofs (sub-control-volume-face based) forms a
 * *closed* subsystem: interior scvf fluxes appear with opposite sign on
 * their two adjacent control volumes and telescope to exactly the boundary
 * flux, independent of what the non-CV dofs are doing. So it is
 * `sum_{vertices} value*Extrusion::volume` -- not the full-domain field
 * integral -- that is the quantity a no-flux (or globally flux-balanced)
 * boundary condition exactly conserves (up to the nonlinear solver's own
 * convergence tolerance). Non-CV dofs' own flux (a Galerkin weak form) does
 * *not* independently telescope to zero the same way, so including them in
 * the sum reintroduces exactly the "leakage" this function avoids.
 *
 * \param gridGeometry The grid geometry
 * \param gridVariablesCache The (current or previous) grid variables cache
 * \param sol The solution vector matching `gridVariablesCache`
 * \param valueAtDof A callable `(vars) -> Scalar` evaluating the integrand
 *        from a control-volume dof's own (nodal) variables.
 */
template<class GridGeometry, class GridVariablesCache, class SolutionVector, class ValueAtDof>
auto integrateControlVolumeField(const GridGeometry& gridGeometry,
                                 const GridVariablesCache& gridVariablesCache,
                                 const SolutionVector& sol,
                                 const ValueAtDof& valueAtDof)
{
    using Scalar = typename SolutionVector::block_type::value_type;
    using Extrusion = Extrusion_t<GridGeometry>;

    auto fvGeometry = localView(gridGeometry);
    auto elemVars = localView(gridVariablesCache);

    Scalar integral = 0.0;
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        fvGeometry.bind(element);
        elemVars.bind(element, fvGeometry, sol);

        for (const auto& scv : scvs(fvGeometry))
            integral += valueAtDof(elemVars[scv])*Extrusion::volume(fvGeometry, scv)*elemVars[scv].extrusionFactor();
    }

    return integral;
}

} // end namespace Dumux::PhaseField

#endif
