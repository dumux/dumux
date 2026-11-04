// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Helper to build Dirichlet constraint lists for CVFE discretizations.
 *
 * Works with problems providing the new boundary condition interface:
 *   `problem.boundaryTypes(fvGeometry, boundaryFace)` returning
 *   `Experimental::BoundaryTypes<numEq>`.
 * Dirichlet DOFs are those for which `!bcTypes.isFluxBoundary(eqIdx)` for some eqIdx.
 */
#ifndef DUMUX_CVFE_APPEND_DIRICHLET_CONSTRAINTS__HH
#define DUMUX_CVFE_APPEND_DIRICHLET_CONSTRAINTS__HH

#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::CVFE {

/*!
 * \brief Append Dirichlet constraints to a constraint list for CVFE discretizations.
 *
 * \param problem      Problem providing `gridGeometry()` and
 *                     `boundaryTypes(fvGeometry, boundaryFace)`.
 * \param getDirichlet Callable `(fvGeometry, boundaryFace, localDof) -> Values`
 *                     returning the Dirichlet values at the DOF.
 * \param constraints  Output vector to append DirichletConstraintData entries to.
 */
template<class Problem, class GetDirichletValues, class ConstraintVector>
void appendDirichletConstraints(const Problem& problem,
                                GetDirichletValues&& getDirichlet,
                                ConstraintVector& constraints)
{
    using ConstraintData = typename ConstraintVector::value_type;
    using ConstraintInfo = typename ConstraintData::ConstraintInfo;
    static constexpr int numEq = ConstraintInfo::size();

    const auto& gg = problem.gridDiscretization();
    auto fvGeometry = localView(gg);

    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bind(element);

        for (const auto& bf : boundaryFaces(fvGeometry))
        {
            const auto bcTypes = problem.boundaryTypes(fvGeometry, bf);

            // Skip if all equations are flux (pure Neumann face)
            if (bcTypes.hasOnlyFluxBoundary())
                continue;

            for (const auto& localDof : localDofs(fvGeometry, bf))
            {
                ConstraintInfo info;
                bool anyDirichlet = false;
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (!bcTypes.isFluxBoundary(static_cast<unsigned>(eqIdx)))
                    {
                        info.set(eqIdx);
                        anyDirichlet = true;
                    }
                }

                if (!anyDirichlet)
                    continue;

                auto values = getDirichlet(fvGeometry, bf, localDof);
                constraints.push_back(ConstraintData{std::move(info), std::move(values), localDof.dofIndex()});
            }
        }
    }
}

} // end namespace Dumux::CVFE
#endif
