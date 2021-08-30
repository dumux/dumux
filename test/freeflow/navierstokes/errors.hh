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
 * \copydoc Dumux::NavierStokesTestError
 */
#ifndef DUMUX_TEST_ERRORS_HH
#define DUMUX_TEST_ERRORS_HH

#include <vector>
#include <cmath>
#include <dumux/discretization/extrusion.hh>
#include <dumux/io/format.hh>

namespace Dumux {
template<class Problem>
class NavierStokesErrorCSVWriter
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Indices = typename Problem::Indices;

    static constexpr int dim = GridGeometry::GridView::dimension;

public:
    NavierStokesErrorCSVWriter()
    { }

    void init(std::string name, int numCCDofs, int numFaceDofs)
    {
        name_ = name;

        //prepare headlines in the error file
        std::ofstream logFile(name_ + "_error.csv");

        logFile << "cc dofs, face dofs, all dofs" << std::endl;

        logFile << Fmt::format("{},{},{}",numCCDofs,numFaceDofs,numCCDofs + numFaceDofs) << std::endl << std::endl;

        logFile << ",absoute L2 error,";
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << ",";
        logFile << "relative L2 error,";
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << ",";
        logFile << "L infinity error,";
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << ",";
        logFile << std::endl;

        logFile << "time";
        for (unsigned int i = 0; i < 3/*absolute L2 error, relative L2error, L infinity error*/; ++i)
        {
            auto uvw = [&](unsigned int dirIdx)
            {
                if (dirIdx == 0)
                    return "u";
                else if (dirIdx == 1)
                    return "v";
                else if (dirIdx == 2)
                    return "w";
                else
                    DUNE_THROW(Dune::InvalidStateException, "dim > 3 not allowed.");
            };

            logFile << ",p";
            for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
                logFile << Fmt::format(",{}",uvw(dirIdx));
        }
        logFile << std::endl;
    }

    template<class PrimaryVariables, class Scalar>
    void printErrors(const PrimaryVariables& l2NormAbs,
                     const PrimaryVariables& l2NormRel,
                     const PrimaryVariables& lInfinityNorm,
                     Scalar time) const
    {
        std::ofstream logFile(name_ + "_error.csv", std::ios::app);
        logFile << Fmt::format("{}", time);

        logFile << Fmt::format(",{}",l2NormAbs[Indices::pressureIdx]);
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::format(",{}",l2NormAbs[Indices::velocity(dirIdx)]);

        logFile << Fmt::format(",{}",l2NormRel[Indices::pressureIdx]);
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::format(",{}",l2NormRel[Indices::velocity(dirIdx)]);

        logFile << Fmt::format(",{}",lInfinityNorm[Indices::pressureIdx]);
        for (unsigned int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::format(",{}",lInfinityNorm[Indices::velocity(dirIdx)]);

        logFile << std::endl;
    }

private:
    std::string name_;
};

template<class Problem, class Scalar = double>
class NavierStokesErrors
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Indices = typename Problem::Indices;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using PrimaryVariables = std::decay_t<decltype(std::declval<Problem>().dirichlet(std::declval<Element>(), std::declval<SubControlVolume>()))>;

    static constexpr int dim = GridGeometry::GridView::dimension;

public:
    NavierStokesErrors(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    {
        std::string name = getParam<std::string>("Problem.Name");
        const int numCCDofs = problem_->gridGeometry().numCellCenterDofs();
        const int numFaceDofs = problem_->gridGeometry().numFaceDofs();
        errorCSVWriter_.init(name, numCCDofs, numFaceDofs);
    }

    template<class SolutionVector>
    void printErrors(const SolutionVector& curSol, Scalar time = 0.0) const
    {
        PrimaryVariables l2NormAbs(0.0);
        PrimaryVariables l2NormRel(0.0);
        PrimaryVariables lInfinityNorm(0.0);

        calculateErrors_(l2NormAbs, l2NormRel, lInfinityNorm, curSol, time);
        errorCSVWriter_.printErrors(l2NormAbs, l2NormRel, lInfinityNorm, time);
    }

private:
    template<class T>
    T absDiff_(const T& a, const T& b) const
    {
        return std::abs(a-b);
    }

    template<class SolutionVector>
    void processPressure_(const SolutionVector& curSol,
                          const SubControlVolume& scv,
                          Scalar& sumPReference,
                          Scalar& sumPError,
                          Scalar& maxPError,
                          Scalar time) const
    {
        const auto analyticalSolutionCellCenter = problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
        const auto numericalSolutionCellCenter = curSol[GridGeometry::cellCenterIdx()][scv.dofIndex()][Indices::pressureIdx - dim];

        const Scalar pError = absDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter);
        const Scalar pReference = absDiff_(analyticalSolutionCellCenter, 0.0);

        maxPError = std::max(maxPError, pError);
        sumPError += pError * pError * Extrusion::volume(scv);
        sumPReference += pReference * pReference * Extrusion::volume(scv);
    }

    template<class SolutionVector>
    void processVelocity_(const SolutionVector& curSol,
                          const SubControlVolume& scv,
                          const SubControlVolumeFace& scvf,
                          Scalar& sumVelReference,
                          Scalar& sumVelError,
                          Scalar& maxVelError,
                          Scalar time) const
    {
        auto faceScvCenter = scv.center() + scvf.center();
        faceScvCenter *= 0.5;
        typename GridGeometry::Traits::FaceSubControlVolume faceScv(faceScvCenter, 0.5*Extrusion::volume(scv));
        const Scalar staggeredHalfVolume = Extrusion::volume(faceScv);

        const auto analyticalSolutionFace = problem_->analyticalSolution(scvf.center(), time)[Indices::velocity(scvf.directionIndex())];
        const auto numericalSolutionFace = curSol[GridGeometry::faceIdx()][scvf.dofIndex()][0];

        const Scalar velError = absDiff_(analyticalSolutionFace, numericalSolutionFace);
        const Scalar velReference = absDiff_(analyticalSolutionFace, 0.0);

        maxVelError = std::max(maxVelError, velError);
        sumVelError += velError * velError * staggeredHalfVolume;
        sumVelReference += velReference * velReference * staggeredHalfVolume;
    }


    /*!
      * \brief Calculates the discrete L2 and Linfinity error between the analytical solution and the numerical approximation.
      *
      * \param problem The object specifying the problem which ought to be simulated
      * \param curSol Vector containing the current solution
      */
    template<class SolutionVector>
    void calculateErrors_(PrimaryVariables& l2NormAbs,
                          PrimaryVariables& l2NormRel,
                          PrimaryVariables& lInfinityNorm,
                          const SolutionVector& curSol,
                          Scalar time) const
    {
        // calculate helping variables
        Scalar totalVolume = 0.;

        PrimaryVariables sumReference(0.0);
        PrimaryVariables sumError(0.0);
        PrimaryVariables maxError(0.0);

        auto fvGeometry = localView(problem_->gridGeometry());
        for (const auto& element : elements(problem_->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                totalVolume += Extrusion::volume(scv);

                processPressure_(curSol, scv, sumReference[Indices::pressureIdx], sumError[Indices::pressureIdx], maxError[Indices::pressureIdx], time);

                // treat face dofs
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    unsigned int index = Indices::velocity(scvf.directionIndex());
                    processVelocity_(curSol, scv, scvf, sumReference[index], sumError[index], maxError[index], time);
                }
            }
        }

        // calculate errors
        for (unsigned int i = 0; i < l2NormAbs.size(); ++i)
        {
            l2NormAbs[i] = std::sqrt(sumError[i] / totalVolume);
            l2NormRel[i] = std::sqrt(sumError[i] / sumReference[i]);
            lInfinityNorm[i] = maxError[i];
        }
    }

    std::shared_ptr<const Problem> problem_;
    NavierStokesErrorCSVWriter<Problem> errorCSVWriter_;
};

template<class Problem>
NavierStokesErrors(std::shared_ptr<Problem> p)
-> NavierStokesErrors<Problem>;

template<class Problem, class Scalar>
NavierStokesErrors(std::shared_ptr<Problem> p, Scalar t)
-> NavierStokesErrors<Problem, Scalar>;
} // end namespace Dumux

#endif
