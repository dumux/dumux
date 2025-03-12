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

#include <unordered_map>

#include <dumux/common/indextraits.hh>

namespace Dumux {

template<class ConstraintInfo, class ConstraintValues>
struct ConstraintData
{
    const auto& constraintInfo() const
    { return info_; }

    const auto& values() const
    {return values_; }

    ConstraintInfo info_;
    ConstraintValues values_;
};

/*!
 * \ingroup Discretization
 * \brief Constraints related to Dirichlet boundaries.
 *        This class generates global constraints based on Dirichlet conditions
 *        set on grid intersections
 */
template<class GG, class CInfo, class CValues>
class DirichletConstraints
{
private:
    using GridIndexType = typename IndexTraits<typename GG::GridView>::GridIndex;
    using ConstraintInfo = CInfo;
    using ConstraintValues = CValues;
public:

    /*!
     * \brief Update the boundary types for all element intersections.
     *
     * \param problem The problem object which needs to be simulated
     * \param dirichletFunction Function for evaluating Dirichlet boundary conditions
     */
    template<class Problem, typename DirichletFunction>
    void update(const Problem& problem, const DirichletFunction& dirichletFunction)
    {
        constraintInfoMap_.clear();

        auto fvGeometry = localView(problem.gridGeometry());
        for (const auto& element : elements(problem.gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& intersection : intersections(problem.gridGeometry().gridView(), element))
            {
                if(!intersection.boundary() || intersection.neighbor())
                    continue;

                const auto& bcTypes = problem.boundaryTypes(fvGeometry, intersection);
                if (bcTypes.hasDirichlet())
                {
                    for (const auto& localDof : localDofs(fvGeometry, intersection))
                    {
                        ConstraintInfo info;
                        // set the Dirichlet constraints
                        for (int eqIdx = 0; eqIdx < ConstraintValues::size(); ++eqIdx)
                            if (bcTypes.isDirichlet(eqIdx))
                                info.set(bcTypes.eqToDirichletIndex(eqIdx), eqIdx);

                        auto dirichletValues = dirichletFunction(fvGeometry, intersection, localDof);
                        constraintInfoMap_[localDof.dofIndex()] = ConstraintData{std::move(info), std::move(dirichletValues)};
                    }
                }
            }
        }
    }

    const auto& map() const
    { return constraintInfoMap_; }

private:
    std::unordered_map<GridIndexType, ConstraintData<ConstraintInfo, ConstraintValues>> constraintInfoMap_;
};


} // namespace Dumux

#endif
