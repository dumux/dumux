// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Constraints related to Dirichlet boundaries
 */
#ifndef DUMUX_DIRICHLET_CONSTRAINTS_HH
#define DUMUX_DIRICHLET_CONSTRAINTS_HH

#include <concepts>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/indextraits.hh>

namespace Dumux {

template<class Info, class Values>
struct ConstraintData
{
    using ConstraintInfo = Info;
    using ConstraintValues = Values;

    const ConstraintInfo& constraintInfo() const
    { return info_; }

    const ConstraintValues& values() const
    { return values_; }

    ConstraintInfo info_;
    ConstraintValues values_;
};

template<class DirichletConstraintInfo, class DirichletValues, class IndexType>
struct DirichletConstraintData : public ConstraintData<DirichletConstraintInfo, DirichletValues>
{
    using GridIndexType = IndexType;

    GridIndexType dofIndex() const
    { return dofIdx_; }

    GridIndexType dofIdx_;
};

/*!
 * \ingroup Discretization
 * \brief Append the Dirichlet constraints of the local dofs on the boundary faces of a grid discretization
 *
 * Each local dof of each boundary face is passed to `constraintInfo(elemDisc, boundaryFace, localDof)`,
 * which returns the constraint info with the equations constrained there (none for a free dof). A dof on
 * several boundary faces gets a single constraint of all equations constrained on any of them. Its values
 * `values(elemDisc, localDof)` are evaluated once, for the first boundary face that constrains the dof.
 * The constraints already in `constraints` are kept unchanged.
 *
 * \param constraints the constraints (e.g. DirichletConstraintData), extended in place
 * \param gridDiscretization the grid discretization
 * \param constraintInfo returns the constraint info of a local dof on a boundary face
 * \param values returns the constraint values of a local dof
 */
template<class Constraint, class GridDiscretization, class ConstraintInfoFunc, class ValuesFunc>
void appendBoundaryDirichletConstraints(std::vector<Constraint>& constraints,
                                        const GridDiscretization& gridDiscretization,
                                        const ConstraintInfoFunc& constraintInfo,
                                        const ValuesFunc& values)
{
    using ConstraintInfo = typename Constraint::ConstraintInfo;
    using ConstraintValues = typename Constraint::ConstraintValues;
    using GridIndexType = typename Constraint::GridIndexType;
    static constexpr auto noConstraint = std::numeric_limits<std::size_t>::max();

    std::vector<std::size_t> constraintIndex(gridDiscretization.numDofs(), noConstraint);

    auto elemDisc = localView(gridDiscretization);
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bindElement(element);
        for (const auto& boundaryFace : boundaryFaces(elemDisc))
        {
            for (const auto& localDof : localDofs(elemDisc, boundaryFace))
            {
                const ConstraintInfo info = constraintInfo(elemDisc, boundaryFace, localDof);
                auto& index = constraintIndex[localDof.dofIndex()];
                for (int eqIdx = 0; eqIdx < ConstraintInfo::size(); ++eqIdx)
                {
                    if (!info.isConstraintEquation(eqIdx))
                        continue;

                    if (index == noConstraint)
                    {
                        index = constraints.size();
                        constraints.push_back(Constraint{ConstraintInfo{},
                                                         ConstraintValues(values(elemDisc, localDof)),
                                                         static_cast<GridIndexType>(localDof.dofIndex())});
                    }

                    auto& mergedInfo = constraints[index].info_;
                    if (!mergedInfo.isConstraintEquation(eqIdx))
                        mergedInfo.set(info.eqToPriVarIndex(eqIdx), eqIdx);
                }
            }
        }
    }
}

} // end namespace Dumux

#endif
