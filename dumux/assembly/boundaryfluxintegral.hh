// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Free functions for adding boundary flux integral contributions to a residual vector.
 */
#ifndef DUMUX_ASSEMBLY_BOUNDARY_FLUX_INTEGRAL_HH
#define DUMUX_ASSEMBLY_BOUNDARY_FLUX_INTEGRAL_HH

#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/concepts.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief Adds boundary flux contributions to the residual for a single FV sub-control volume face.
 *
 * \param residual The element residual vector
 * \param problem The problem providing boundaryFlux(...)
 * \param elemDisc The bound element discretization (local view)
 * \param elemVars The element variables
 * \param scvf The boundary sub-control volume face
 * \param bcTypes The boundary condition types for this face
 */
template<class ResidualVector, class Problem, Concepts::FVElementDiscretization ElementDiscretization,
         class ElementVariables, class BoundaryTypes>
void addFVBoundaryFluxIntegral(ResidualVector& residual,
                               const Problem& problem,
                               const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const typename ElementDiscretization::SubControlVolumeFace& scvf,
                               const BoundaryTypes& bcTypes)
{
    using BoundaryFluxes = typename Problem::Traits::NumEqVector;
    BoundaryFluxes boundaryFlux(0.0);
    for (const auto& qpData : Dumux::CVFE::quadratureRule(elemDisc, scvf))
        boundaryFlux += qpData.weight() * problem.boundaryFlux(elemDisc, elemVars, qpData.ipData());

    const auto& insideScv = elemDisc.scv(scvf.insideScvIdx());
    boundaryFlux *= elemVars[insideScv].extrusionFactor();

    for (int eqIdx = 0; eqIdx < BoundaryFluxes::dimension; ++eqIdx)
        if (bcTypes.isFluxBoundary(eqIdx))
            residual[insideScv.localDofIndex()][eqIdx] += boundaryFlux[eqIdx];
}

/*!
 * \ingroup Assembly
 * \brief Adds boundary flux contributions to the residual related to FE discretization.
 *
 * \param residual The element residual vector
 * \param problem The problem providing boundaryFlux(...)
 * \param elemDisc The bound element discretization (local view)
 * \param elemVars The element variables
 * \param boundaryFace The boundary face
 * \param bcTypes The boundary condition types for this face
 */
template<class ResidualVector, class Problem, class ElementDiscretization,
         class ElementVariables, class BoundaryTypes>
void addFEBoundaryFluxIntegral(ResidualVector& residual,
                               const Problem& problem,
                               const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const typename ElementDiscretization::BoundaryFace& boundaryFace,
                               const BoundaryTypes& bcTypes)
{
    using BoundaryFluxes = typename Problem::Traits::NumEqVector;

    for (const auto& qpData : Dumux::CVFE::quadratureRule(elemDisc, boundaryFace))
    {
        const auto& ipData = qpData.ipData();
        const auto& ipCache = cache(elemVars, ipData);
        const BoundaryFluxes boundaryFlux = qpData.weight() * problem.boundaryFlux(elemDisc, elemVars, ipData);
        const auto& shapeValues = ipCache.shapeValues();

        for (const auto& localDof : nonCVLocalDofs(elemDisc))
            for (int eqIdx = 0; eqIdx < BoundaryFluxes::dimension; ++eqIdx)
                if (bcTypes.isFluxBoundary(eqIdx))
                    residual[localDof.index()][eqIdx] += shapeValues[localDof.index()] * boundaryFlux[eqIdx];
    }
}

} // namespace Dumux::Experimental

#endif // DUMUX_ASSEMBLY_BOUNDARY_FLUX_INTEGRAL_HH
