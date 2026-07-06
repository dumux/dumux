// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief Shared finite-element residual assembly helpers for local dofs
 *        without an associated sub-control volume (e.g. the edge dofs of the
 *        PQ2/PQ3 hybrid CVFE schemes), reused across phase-field models.
 */
#ifndef DUMUX_PHASEFIELD_COMMON_FE_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_COMMON_FE_LOCAL_RESIDUAL_HH

#include <type_traits>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>

namespace Dumux::PhaseField {

/*!
 * \ingroup PhaseFieldModels
 * \brief Computes \f$ \int N_i(x)\,dx \f$ for every local dof i of the
 *        element (control-volume dofs and non-control-volume dofs alike).
 *
 * These integrals are the diagonal entries of a *consistent* (row-sum)
 * lumped mass matrix: since Lagrange shape functions form a partition of
 * unity (\f$ \sum_i N_i(x) = 1 \f$ everywhere), they always sum to exactly
 * the element's volume, \f$ \sum_i \int N_i\,dx = \int 1\,dx \f$, for any
 * element order -- unlike control-volume dofs' classic sub-control-volume
 * geometric volume (see e.g. `Extrusion::volume(fvGeometry, scv)`), which is
 * a *different* quantity that does not, in general, equal
 * \f$ \int N_i\,dx \f$ for higher-order (e.g. PQ2) bases.
 *
 * \warning Despite that partition-of-unity property, do *not* substitute
 * this weight for control-volume dofs' own storage term (in place of
 * `Extrusion::volume`): empirically, doing so pairs a much smaller
 * "capacitance" at vertex dofs with the (unchanged) much larger sub-control-
 * volume region their flux is collected over, which showed up as first-
 * Newton-iteration relative shifts of several hundred percent (vs. tens of
 * percent normally) and *worse*, not better, mass-conservation/energy-
 * dissipation behavior in practice for the Cahn-Hilliard test. A correct
 * fix would need the flux discretization at control-volume dofs to become
 * consistent with the smaller effective region too (e.g. by redefining the
 * PQ2 hybrid scheme's sub-control-volume geometry itself, in
 * `dumux/discretization/pq2/geometryhelper.hh`, which currently reuses the
 * classic Box scheme's geometry unmodified for corner scvs) -- a shared-
 * framework change, not a local one. This function is only used for the
 * finite-element (non-CV) local dofs' own storage term, below.
 *
 * \param fvGeometry The finite-volume geometry of the element
 */
template<class FVElementGeometry>
auto computeShapeFunctionIntegrals(const FVElementGeometry& fvGeometry)
{
    const auto& localBasis = fvGeometry.feLocalBasis();
    using LocalBasis = std::decay_t<decltype(localBasis)>;
    using RangeType = typename LocalBasis::Traits::RangeType;

    ReservedBlockVector<RangeType, Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>()>
        integralShapeFunctions(localBasis.size());

    const auto& geometry = fvGeometry.elementGeometry();
    const auto& element = fvGeometry.element();
    using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
    using FeIpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;

    for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
    {
        const auto& ipData = qpData.ipData();
        FeIpData feIpData(geometry, ipData.local(), ipData.global(), localBasis);

        for (const auto& localDof : localDofs(fvGeometry))
            integralShapeFunctions[localDof.index()] += qpData.weight() * feIpData.shapeValue(localDof.index());
    }

    return integralShapeFunctions;
}

/*!
 * \ingroup PhaseFieldModels
 * \brief Adds a mass-lumped storage (time-derivative) residual contribution
 *        for local dofs without an associated sub-control volume.
 *
 * Mass-lumping is applied by precomputing the integral of each non-CV local
 * dof's shape function over the element (see computeShapeFunctionIntegrals()),
 * then multiplying by the (possibly multi-equation) time derivative of the
 * quantity returned by \a storageAtDof. A no-op for discretizations without
 * non-CV local dofs (e.g. Box/PQ1Bubble).
 *
 * \param residual The element residual vector to add to
 * \param fvGeometry The finite-volume geometry of the element
 * \param prevElemVars The variables for all local dofs of the element at the previous time level
 * \param curElemVars The variables for all local dofs of the element at the current time level
 * \param timeStepSize The current time step size
 * \param storageAtDof A function taking a local dof's variables and returning
 *                      its (per-equation) storage term, e.g.
 *                      `[](const auto& vars) { NumEqVector s(0.0); s[eqIdx] = vars.someQuantity(); return s; }`
 */
template<class ResidualVector, class FVElementGeometry, class ElementVariables, class Scalar, class StorageAtDof>
void addMassLumpedStorageResidual(ResidualVector& residual,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVariables& prevElemVars,
                                  const ElementVariables& curElemVars,
                                  const Scalar timeStepSize,
                                  StorageAtDof&& storageAtDof)
{
    if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
    {
        // Make sure we don't iterate over quadrature points if there are no non-CV dofs
        if (nonCVLocalDofs(fvGeometry).empty())
            return;

        const auto integralShapeFunctions = computeShapeFunctionIntegrals(fvGeometry);

        for (const auto& localDof : nonCVLocalDofs(fvGeometry))
        {
            const auto localDofIdx = localDof.index();
            const auto storageCur = storageAtDof(curElemVars[localDofIdx]);
            const auto storagePrev = storageAtDof(prevElemVars[localDofIdx]);
            const auto timeDeriv = (storageCur - storagePrev)/timeStepSize;
            residual[localDofIdx] += integralShapeFunctions[localDofIdx][0]*timeDeriv;
        }
    }
}

} // end namespace Dumux::PhaseField

#endif
