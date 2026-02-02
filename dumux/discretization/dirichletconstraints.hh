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
#include <unordered_map>

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

} // end namespace Dumux


namespace Dumux::CVFE::Detail {

template <typename C>
concept DirichletConstraintContainer = requires(
    C& c,
    typename C::value_type v,
    typename C::value_type::ConstraintInfo info,
    typename C::value_type::ConstraintValues vals,
    typename C::value_type::GridIndexType idx
){
    c.push_back(v);
    { typename C::value_type{info, vals, idx} };
};

} // end namespace Dumux::CVFE::Detail


namespace Dumux::CVFE {

/*!
 * \brief Append constraints for Dirichlet boundaries.
 *
 * \param problem The problem object which needs to be simulated
 * \param dirichletFunction Function for evaluating Dirichlet boundary conditions
 * \param constraints The constraints to be appended to
 */
template<class Problem, class DirichletFunction, Detail::DirichletConstraintContainer Constraints>
void appendDirichletConstraints(const Problem& problem,
                                const DirichletFunction& dirichletFunction,
                                Constraints& constraints)
{
    using Data = typename Constraints::value_type;
    using ConstraintInfo = typename Data::ConstraintInfo;
    using DirichletValues = typename Data::ConstraintValues;

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
                    for (int eqIdx = 0; eqIdx < DirichletValues::size(); ++eqIdx)
                        if (bcTypes.isDirichlet(eqIdx))
                            info.set(bcTypes.eqToDirichletIndex(eqIdx), eqIdx);

                    auto dirichletValues = dirichletFunction(fvGeometry, intersection, localDof);
                    constraints.push_back(Data{std::move(info), std::move(dirichletValues), localDof.dofIndex()});
                }
            }
        }
    }
}

} // end namespace Dumux::CVFE

#endif
