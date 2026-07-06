// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawModel
 * \brief Helper functions for assembling the finite-element part of the
 *        hybrid CV/FE Hele-Shaw local residual (contributions of local dofs
 *        without an associated sub-control volume, e.g. the edge dofs of PQ2/PQ3).
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_FE_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_FE_LOCAL_RESIDUAL_HH

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/phasefield/common/felocalresidual.hh>

#include "flux.hh"

namespace Dumux {

/*!
 * \ingroup HeleShawModel
 * \brief Helper class for evaluating the finite-element part of the hybrid
 *        CVFE Hele-Shaw local residual (non-control-volume local dofs).
 *
 * The advective (Darcy) term here is a plain Galerkin (non-upwinded)
 * discretization: upwinding is inherently a face-based concept, and
 * non-control-volume local dofs (e.g. PQ2's edge dofs) have no associated
 * sub-control-volume face to upwind across.
 */
template<class NumEqVector, class Indices>
class HeleShawTwoPCVFEFELocalResidualTerms
{
public:
    /*!
     * \brief Add storage residual contributions for non-CV local dofs.
     *        Only the phase field is transported.
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
                const HeleShawTwoPCVFEFluxFunctionContext context(fvGeometry, elemVars, ipCache);
                const auto& shapeValues = ipCache.shapeValues();

                const auto phi = context.phaseField();
                const auto mu = context.chemicalPotential();
                const auto rhoMix = problem.mixtureDensity(phi);
                const auto etaMix = problem.mixtureViscosity(phi);
                const auto lambda = problem.permeability()/etaMix;
                const auto& g = problem.gravity();

                // Darcy velocity vector v = -lambda*(gradP - mu*gradPhi - rhoMix*g)
                auto velocity = context.gradPressure();
                velocity.axpy(-mu, context.gradPhaseField());
                velocity.axpy(-rhoMix, g);
                velocity *= -lambda;

                auto sourceAtIp = problem.source(fvGeometry, elemVars, ipData);
                // model-intrinsic algebraic identity contribution, closing
                // mu = df/dphi - gamma laplace(phi) (see HeleShawTwoPCVFELocalResidual::sourceIntegral)
                sourceAtIp[Indices::chemPotEqIdx] += mu;

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    const auto& gradN = ipCache.gradN(localDofIdx);
                    NumEqVector fluxAndSourceTerm(0.0);

                    // Darcy continuity (weak form of div(v)=0)
                    fluxAndSourceTerm[Indices::continuityEqIdx] += -(velocity*gradN);

                    // phase field: Galerkin (non-upwinded) advection + Cahn-Hilliard diffusion
                    fluxAndSourceTerm[Indices::phaseFieldEqIdx] +=
                        -phi*(velocity*gradN) + problem.chMobility()*(context.gradChemicalPotential()*gradN);

                    // chemical potential: Laplacian term of the Cahn-Hilliard equation
                    fluxAndSourceTerm[Indices::chemPotEqIdx] += problem.surfaceTensionCoeff()*(context.gradPhaseField()*gradN);

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
