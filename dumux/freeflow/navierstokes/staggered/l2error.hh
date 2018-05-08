// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_L2_ERROR_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_L2_ERROR_HH

#include <vector>
#include <cmath>

// #include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesTests
 * \brief  Routines to calculate the discrete L2 error
 */
template<class Scalar, class ModelTraits, class PrimaryVariables, int dim, bool verbose = false>
class NavierStokesTestL2Error
{
    using Indices = typename ModelTraits::Indices;

public:
    /*!
      * \brief Calculate the L2 error between the analytical solution and the numerical approximation.
      *
      * \param problem The problem
      * \param curSol Vector containing the current solution
      * \param pOrder Polynomial order of quadrature rule
      */
    template<class Problem, class SolutionVector>
    static void printL2Error(const Problem& problem, const SolutionVector& curSol, int pOrder = 1)
    {
        const auto l2error = calculateL2Error(problem, curSol, pOrder);
        std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                  << std::setw(6) << problem.fvGridGeometry().numCellCenterDofs() << " cc dofs and "
                  << problem.fvGridGeometry().numFaceDofs() << " face dofs (total: "
                  << problem.fvGridGeometry().numCellCenterDofs() + problem.fvGridGeometry().numFaceDofs() << "): "
                  << std::scientific
                  << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx];
        if (verbose)
        {
            std::cout << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                      << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx];
        }
        else
        {

            std::cout << " , L2(v) = " << l2error.first[Indices::velocity(0)] << " / " << l2error.second[Indices::velocity(0)];
        }
        std::cout << std::endl;
    }

    /*!
      * \brief Calculate the L2 error between the analytical solution and the numerical approximation.
      *
      * \param problem The problem
      * \param curSol Vector containing the current solution
      * \param pOrder Polynomial order of quadrature rule
      */
    template<class Problem, class SolutionVector>
    static auto calculateL2Error(const Problem& problem, const SolutionVector& curSol, int pOrder = 1)
    {
        using FVGridGeometry = std::decay_t<decltype(problem.fvGridGeometry())>;
        PrimaryVariables sumError(0.0), sumReference(0.0), totalVolume(0.0), l2NormAbs(0.0), l2NormRel(0.0);

        const int numFaceDofs = problem.fvGridGeometry().numFaceDofs();

        std::vector<Scalar> staggeredVolume(numFaceDofs);
        std::vector<Scalar> errorVelocity(numFaceDofs);
        std::vector<Scalar> velocityReference(numFaceDofs);
        std::vector<int> directionIndex(numFaceDofs);

        for (const auto& element : elements(problem.fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(problem.fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                // integrate over element using a quadrature rule
                const auto& rule = Dune::QuadratureRules<Scalar, dim>::rule(element.geometry().type(), pOrder);
                for (typename Dune::QuadratureRule<Scalar, dim>::const_iterator qIt = rule.begin(); qIt != rule.end(); ++qIt)
                {
                    auto localPos = qIt->position();
                    const auto dofIdxCellCenter = scv.dofIndex();
                    const auto analyticalSolutionCellCenter = problem.dirichletAtPos(element.geometry().global(localPos))[Indices::pressureIdx];
                    const auto numericalSolutionCellCenter = curSol[FVGridGeometry::cellCenterIdx()][dofIdxCellCenter][Indices::pressureIdx - ModelTraits::dim()];
                    Scalar weigth = qIt->weight() * element.geometry().integrationElement(localPos);
                    sumError[Indices::pressureIdx] += squaredDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter) * weigth;
                    sumReference[Indices::pressureIdx] += analyticalSolutionCellCenter * analyticalSolutionCellCenter * weigth;
                    totalVolume[Indices::pressureIdx] += weigth;

                    // treat face dofs
                    for (auto&& scvf : scvfs(fvGeometry))
                    {
                        const int dofIdxFace = scvf.dofIndex();
                        const int dirIdx = scvf.directionIndex();
                        // only treat the part of the quadrature which belongs to this dof
                        if (((scvf.localFaceIdx() % 2 == 0) && localPos[dirIdx] > 0.5)
                             || ((scvf.localFaceIdx() % 2 == 1) && localPos[dirIdx] < 0.5 + 1e-10))
                            continue;

                        const auto analyticalSolutionFace = problem.dirichletAtPos(element.geometry().global(localPos))[Indices::velocity(dirIdx)];
                        const auto numericalSolutionFace = curSol[FVGridGeometry::faceIdx()][dofIdxFace][0];
                        directionIndex[dofIdxFace] = dirIdx;
                        errorVelocity[dofIdxFace] += squaredDiff_(analyticalSolutionFace, numericalSolutionFace) * weigth;
                        velocityReference[dofIdxFace] += squaredDiff_(analyticalSolutionFace, 0.0) * weigth;
                        staggeredVolume[dofIdxFace] += weigth;
                    }
                }
            }
        }

        // get the absolute and relative discrete L2-error for cell-center dofs
        l2NormAbs[Indices::pressureIdx] = std::sqrt(sumError[Indices::pressureIdx] / totalVolume[Indices::pressureIdx]);
        l2NormRel[Indices::pressureIdx] = std::sqrt(sumError[Indices::pressureIdx] / sumReference[Indices::pressureIdx]);

        // get the absolute and relative discrete L2-error for face dofs
        for(int i = 0; i < numFaceDofs; ++i)
        {
            const auto error = errorVelocity[i];
            const auto ref = velocityReference[i];
            if (verbose)
            {
                sumError[Indices::velocity(directionIndex[i])] += error;
                sumReference[Indices::velocity(directionIndex[i])] += ref;
                totalVolume[Indices::velocity(directionIndex[i])] += staggeredVolume[i];
            }
            else
            {
                sumError[Indices::velocity(0)] += error;
                sumReference[Indices::velocity(0)] += ref;
                totalVolume[Indices::velocity(0)] += staggeredVolume[i];
            }
        }

        for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
        {
            l2NormAbs[Indices::velocity(dirIdx)] = std::sqrt(sumError[Indices::velocity(dirIdx)] / totalVolume[Indices::velocity(dirIdx)]);
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
}

#endif
