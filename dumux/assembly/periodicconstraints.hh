// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Periodic constraints of assembled linear systems
 *
 * Of each set of periodically identified degrees of freedom, the one with the smallest index
 * keeps the sum of the equations of the set, and the equations of the others are replaced by
 * the constraint that their solution equals the one of the smallest index.
 */
#ifndef DUMUX_ASSEMBLY_PERIODIC_CONSTRAINTS_HH
#define DUMUX_ASSEMBLY_PERIODIC_CONSTRAINTS_HH

#include <algorithm>
#include <cassert>

#include <dumux/common/deprecated.hh>
#include <dumux/common/typetraits/periodic.hh>

namespace Dumux::Detail::PeriodicConstraints {

//! Call f(first, second) for each periodically mapped degree of freedom second of the degree of freedom first with the smallest index
template<class GridGeometry, class Function>
void forEachPeriodicPair(const GridGeometry& gridGeometry, Function&& f)
{
    if constexpr (Dumux::Detail::hasPeriodicDofMap<GridGeometry>())
    {
        for (const auto& m : gridGeometry.periodicDofMap())
        {
            const auto& periodicallyMappedDofs = Dumux::Deprecated::ensureRangeOfPeriodicDofs(m.second);
            if (std::ranges::all_of(periodicallyMappedDofs, [&](auto second){ return m.first < second; }))
                for (const auto& second : periodicallyMappedDofs)
                    f(m.first, second);
        }
    }
}

//! Add the residual of each periodic degree of freedom to the one it is mapped to and replace it by the constraint
template<class GridGeometry, class Residual, class Solution>
void enforceInResidual(const GridGeometry& gridGeometry, Residual& res, const Solution& curSol)
{
    forEachPeriodicPair(gridGeometry, [&](const auto first, const auto second)
    {
        res[first] += res[second];
        res[second] = curSol[second] - curSol[first];
    });
}

//! Add the matrix row of each periodic degree of freedom to the one it is mapped to and replace it by the derivative of the constraint
template<class GridGeometry, class Matrix>
void enforceInJacobian(const GridGeometry& gridGeometry, Matrix& jac)
{
    forEachPeriodicPair(gridGeometry, [&](const auto first, const auto second)
    {
        const auto end = jac[second].end();
        for (auto it = jac[second].begin(); it != end; ++it)
            jac[first][it.index()] += (*it);

        for (auto it = jac[second].begin(); it != end; ++it)
        {
            auto& matrixBlock = *it;
            matrixBlock = 0.0;

            assert(matrixBlock.N() == matrixBlock.M());
            const double diagValue = it.index() == second ? 1.0 : (it.index() == first ? -1.0 : 0.0);
            if (diagValue != 0.0)
                for (int eIdx = 0; eIdx < matrixBlock.N(); ++eIdx)
                    matrixBlock[eIdx][eIdx] = diagValue;
        }
    });
}

/*!
 * \brief Add the row of each periodic degree of freedom to the one it is mapped to in a coupling block
 *        of a multidomain matrix (the rows of the periodic domain) and zero it, as the constraint
 *        does not depend on other domains
 */
template<class GridGeometry, class Matrix>
void enforceInCouplingJacobian(const GridGeometry& gridGeometry, Matrix& jacCoupling)
{
    forEachPeriodicPair(gridGeometry, [&](const auto first, const auto second)
    {
        const auto end = jacCoupling[second].end();
        for (auto it = jacCoupling[second].begin(); it != end; ++it)
            jacCoupling[first][it.index()] += (*it);

        for (auto it = jacCoupling[second].begin(); it != end; ++it)
            (*it) = 0.0;
    });
}

} // end namespace Dumux::Detail::PeriodicConstraints

#endif
