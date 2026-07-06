// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldLevelSetModel
 * \brief Helper functions for assembling the finite-element part of the
 *        hybrid CV/FE local residual (contributions of local dofs without
 *        an associated sub-control volume, e.g. the edge dofs of PQ2/PQ3).
 */
#ifndef DUMUX_PHASEFIELD_LEVEL_SET_FE_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_LEVEL_SET_FE_LOCAL_RESIDUAL_HH

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/phasefield/common/felocalresidual.hh>

#include "flux.hh"

namespace Dumux {

/*!
 * \ingroup PhaseFieldLevelSetModel
 * \brief Helper class for evaluating the finite-element part of the
 *        conservative phase-field level-set local residual (non-control-volume local dofs).
 */
template<class NumEqVector, class Indices>
class PhaseFieldLevelSetFELocalResidualTerms
{
public:
    /*!
     * \brief Add storage residual contributions for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current time level
     * \param timeStepSize The current time step size
     */
    template<class ResidualVector, class Problem, class FVElementGeometry, class ElementVariables, class Scalar>
    static void addStorageTerms(ResidualVector& residual,
                                const Problem& problem,
                                const FVElementGeometry& fvGeometry,
                                const ElementVariables& prevElemVars,
                                const ElementVariables& curElemVars,
                                const Scalar timeStepSize)
    {
        PhaseField::addMassLumpedStorageResidual(residual, fvGeometry, prevElemVars, curElemVars, timeStepSize,
            [](const auto& vars)
            {
                NumEqVector storage(0.0);
                storage[Indices::phaseFieldEqIdx] = vars.phaseField();
                return storage;
            });
    }

    /*!
     * \brief Add flux and source residual contributions for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     */
    template<class ResidualVector, class Problem, class FVElementGeometry, class ElementVariables>
    static void addFluxAndSourceTerms(ResidualVector& residual,
                                      const Problem& problem,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVariables& elemVars)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            // Make sure we don't iterate over quadrature points if there are no non-CV dofs
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            const auto& element = fvGeometry.element();
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const LevelSetFluxFunctionContext context(fvGeometry, elemVars, ipCache);
                const auto& shapeValues = ipCache.shapeValues();

                const auto phi = context.phaseField();
                const auto compressionCoeff = phi*(1.0 - phi);
                const auto normal = context.normal();

                const auto sourceAtIp = problem.source(fvGeometry, elemVars, ipData);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    const auto& gradN = ipCache.gradN(localDofIdx);
                    NumEqVector fluxAndSourceTerm(0.0);

                    // weak-form terms: -grad(N_i) . (phi(1-phi) n) [compression] + grad(N_i) . (epsilon grad(phi)) [diffusion]
                    fluxAndSourceTerm[Indices::phaseFieldEqIdx] +=
                        -compressionCoeff*(normal*gradN) + problem.epsilon()*(context.gradPhaseField()*gradN);

                    // reaction/source term and scatter to residual
                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    {
                        fluxAndSourceTerm[eqIdx] -= shapeValues[localDofIdx] * sourceAtIp[eqIdx];
                        residual[localDofIdx][eqIdx] += qpData.weight()*fluxAndSourceTerm[eqIdx];
                    }
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
