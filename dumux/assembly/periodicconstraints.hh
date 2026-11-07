// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Enforce periodicity in an assembled linear system
 *
 * The degrees of freedom on periodic boundaries exist once per side. The row of the second
 * dof is added to the row of the first, which yields the balance over the joined control
 * volume or the joined basis function, and is replaced by the constraint that both dofs
 * have the same value.
 */
#ifndef DUMUX_ASSEMBLY_PERIODIC_CONSTRAINTS_HH
#define DUMUX_ASSEMBLY_PERIODIC_CONSTRAINTS_HH

#include <algorithm>
#include <cassert>
#include <ranges>

#include <dumux/common/deprecated.hh>
#include <dumux/common/typetraits/periodic.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Assembly
 * \brief Call f(first, second) for each pair of periodically mapped dofs,
 *        with first being the lowest index of all dofs mapped onto each other
 */
template<class GridDiscretization, class Function>
void forEachPeriodicDofPair(const GridDiscretization& gridDiscretization, Function&& f)
{
    if constexpr (hasPeriodicDofMap<GridDiscretization>())
    {
        for (const auto& m : gridDiscretization.periodicDofMap())
        {
            const auto& periodicallyMappedDofs = Dumux::Deprecated::ensureRangeOfPeriodicDofs(m.second);
            if (std::ranges::all_of(periodicallyMappedDofs, [&](auto second){ return m.first < second; }))
                for (const auto& second : periodicallyMappedDofs)
                    f(m.first, second);
        }
    }
}

/*!
 * \ingroup Assembly
 * \brief Move the rows of the second dofs of periodic pairs to the first dofs in the Jacobian
 *        and write the derivatives of the periodicity constraint into the rows of the second dofs
 */
template<class JacobianMatrix, class GridDiscretization>
void enforcePeriodicConstraints(JacobianMatrix& jac, const GridDiscretization& gridDiscretization)
{
    forEachPeriodicDofPair(gridDiscretization, [&](const auto first, const auto second)
    {
        const auto end = jac[second].end();
        for (auto it = jac[second].begin(); it != end; ++it)
            jac[first][it.index()] += (*it);

        // the constraint has the derivatives id for second and -id for first
        auto setMatrixBlock = [] (auto& matrixBlock, double diagValue)
        {
            for (int eIdx = 0; eIdx < matrixBlock.N(); ++eIdx)
                matrixBlock[eIdx][eIdx] = diagValue;
        };

        for (auto it = jac[second].begin(); it != end; ++it)
        {
            auto& matrixBlock = *it;
            matrixBlock = 0.0;

            assert(matrixBlock.N() == matrixBlock.M());
            if (it.index() == second)
                setMatrixBlock(matrixBlock, 1.0);

            if (it.index() == first)
                setMatrixBlock(matrixBlock, -1.0);
        }
    });
}

/*!
 * \ingroup Assembly
 * \brief Move the entries of the second dofs of periodic pairs to the first dofs in the residual
 *        and write the periodicity constraint into the entries of the second dofs
 */
template<class ResidualType, class SolutionVector, class GridDiscretization>
void enforcePeriodicConstraints(ResidualType& res, const SolutionVector& curSol, const GridDiscretization& gridDiscretization)
{
    forEachPeriodicDofPair(gridDiscretization, [&](const auto first, const auto second)
    {
        res[first] += res[second];
        res[second] = curSol[second] - curSol[first];
    });
}

/*!
 * \ingroup Assembly
 * \brief Enforce the periodicity constraints in Jacobian and residual
 */
template<class JacobianMatrix, class ResidualType, class SolutionVector, class GridDiscretization>
void enforcePeriodicConstraints(JacobianMatrix& jac, ResidualType& res, const SolutionVector& curSol, const GridDiscretization& gridDiscretization)
{
    enforcePeriodicConstraints(jac, gridDiscretization);
    enforcePeriodicConstraints(res, curSol, gridDiscretization);
}

} // end namespace Dumux::Detail

#endif
