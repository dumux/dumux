// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesTests
 * \copydoc Dumux::NavierStokesTestL2Error
 */
#ifndef DUMUX_TEST_L2_ERROR_HH
#define DUMUX_TEST_L2_ERROR_HH

#include <vector>
#include <cmath>

namespace Dumux
{

/*!
 * \ingroup NavierStokesTests
 * \brief  Routines to calculate the discrete L2 error
 */
template<class Scalar, class ModelTraits, class PrimaryVariables>
class NavierStokesTestL2Error
{
    using Indices = typename ModelTraits::Indices;

public:

    /*!
      * \brief Calculate the L2 error between the analytical solution and the numerical approximation.
      *
      * \param problem The object specifying the problem which ought to be simulated
      * \param curSol Vector containing the current solution
      */
    template<class Problem, class SolutionVector>
    static auto calculateL2Error(const Problem& problem, const SolutionVector& curSol)
    {
        using FVGridGeometry = std::decay_t<decltype(problem.fvGridGeometry())>;
        PrimaryVariables sumError(0.0), sumReference(0.0), l2NormAbs(0.0), l2NormRel(0.0);

        const int numFaceDofs = problem.fvGridGeometry().numFaceDofs();

        std::vector<Scalar> staggeredVolume(numFaceDofs);
        std::vector<Scalar> errorVelocity(numFaceDofs);
        std::vector<Scalar> velocityReference(numFaceDofs);
        std::vector<int> directionIndex(numFaceDofs);

        Scalar totalVolume = 0.0;

        for (const auto& element : elements(problem.fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(problem.fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                // treat cell-center dofs
                const auto dofIdxCellCenter = scv.dofIndex();
                const auto& posCellCenter = scv.dofPosition();
                const auto analyticalSolutionCellCenter = problem.dirichletAtPos(posCellCenter)[Indices::pressureIdx];
                const auto numericalSolutionCellCenter = curSol[FVGridGeometry::cellCenterIdx()][dofIdxCellCenter][Indices::pressureIdx - ModelTraits::dim()];
                sumError[Indices::pressureIdx] += squaredDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter) * scv.volume();
                sumReference[Indices::pressureIdx] += analyticalSolutionCellCenter * analyticalSolutionCellCenter * scv.volume();
                totalVolume += scv.volume();

                // treat face dofs
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const int dofIdxFace = scvf.dofIndex();
                    const int dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionFace = problem.dirichletAtPos(scvf.center())[Indices::velocity(dirIdx)];
                    const auto numericalSolutionFace = curSol[FVGridGeometry::faceIdx()][dofIdxFace][0];
                    directionIndex[dofIdxFace] = dirIdx;
                    errorVelocity[dofIdxFace] = squaredDiff_(analyticalSolutionFace, numericalSolutionFace);
                    velocityReference[dofIdxFace] = squaredDiff_(analyticalSolutionFace, 0.0);
                    const Scalar staggeredHalfVolume = 0.5 * scv.volume();
                    staggeredVolume[dofIdxFace] = staggeredVolume[dofIdxFace] + staggeredHalfVolume;
                }
            }
        }

        // get the absolute and relative discrete L2-error for cell-center dofs
        l2NormAbs[Indices::pressureIdx] = std::sqrt(sumError[Indices::pressureIdx] / totalVolume);
        l2NormRel[Indices::pressureIdx] = std::sqrt(sumError[Indices::pressureIdx] / sumReference[Indices::pressureIdx]);

        // get the absolute and relative discrete L2-error for face dofs
        for(int i = 0; i < numFaceDofs; ++i)
        {
            const int dirIdx = directionIndex[i];
            const auto error = errorVelocity[i];
            const auto ref = velocityReference[i];
            const auto volume = staggeredVolume[i];
            sumError[Indices::velocity(dirIdx)] += error * volume;
            sumReference[Indices::velocity(dirIdx)] += ref * volume;
        }

        for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
        {
            l2NormAbs[Indices::velocity(dirIdx)] = std::sqrt(sumError[Indices::velocity(dirIdx)] / totalVolume);
            l2NormRel[Indices::velocity(dirIdx)] = std::sqrt(sumError[Indices::velocity(dirIdx)] / sumReference[Indices::velocity(dirIdx)]);
        }
        return std::make_pair(l2NormAbs, l2NormRel);
    }

private:

    template<class T>
    static T squaredDiff_(const T& a, const T& b)
    {
        return (a-b)*(a-b);
    }

};
} //end namespace DUMUX_TEST_L2_ERROR_HH

#endif
