// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardModel
 * \brief Helper functions for assembling the finite-element part of the
 *        hybrid CV/FE local residual (contributions of local dofs without
 *        an associated sub-control volume, e.g. the edge dofs of PQ2/PQ3).
 */
#ifndef DUMUX_PHASEFIELD_CAHNHILLIARD_FE_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_CAHNHILLIARD_FE_LOCAL_RESIDUAL_HH

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/phasefield/common/felocalresidual.hh>

#include "flux.hh"

namespace Dumux {

/*!
 * \ingroup CahnHilliardModel
 * \brief Helper class for evaluating the finite-element part of the
 *        Cahn-Hilliard local residual (non-control-volume local dofs).
 */
template<class NumEqVector, class Indices>
class CahnHilliardFELocalResidualTerms
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
                // the chemical potential equation has no storage term
                NumEqVector storage(0.0);
                storage[Indices::massBalanceEqIdx] = vars.concentration();
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

            // The double-well reaction term f'(c) is cubic in c (c itself
            // degree-2 on PQ2), so as a function of position it is degree 6;
            // combined with a degree-2 test function this needs an exact
            // degree-8 rule -- the test's properties.hh raises the grid
            // geometry's default element-quadrature order accordingly (the
            // grid variables cache is sized to that scheme, so the order
            // can't be overridden ad-hoc at this call site).
            const auto& element = fvGeometry.element();
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const CahnHilliardFluxFunctionContext context(fvGeometry, elemVars, ipCache);
                const auto& shapeValues = ipCache.shapeValues();

                auto sourceAtIp = problem.source(fvGeometry, elemVars, ipData);
                // model-intrinsic algebraic identity contribution, closing
                // mu = df/dc - gamma laplace(c) (see CahnHilliardLocalResidual::sourceIntegral)
                sourceAtIp[Indices::chemicalPotentialEqIdx] += context.chemicalPotential();

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    NumEqVector fluxAndSourceTerm(0.0);

                    // weak-form diffusive terms: +grad(N_i) . (M grad(mu)) and +grad(N_i) . (gamma grad(c))
                    fluxAndSourceTerm[Indices::massBalanceEqIdx] += problem.mobility()*(context.gradChemicalPotential()*ipCache.gradN(localDofIdx));
                    fluxAndSourceTerm[Indices::chemicalPotentialEqIdx] += problem.surfaceTension()*(context.gradConcentration()*ipCache.gradN(localDofIdx));

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
