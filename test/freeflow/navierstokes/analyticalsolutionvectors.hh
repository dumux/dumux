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
 * \copydoc Dumux::NavierStokesAnalyticalSolutionVectors
 */
#ifndef DUMUX_TEST_ANALYTICALSOLVECTORS_HH
#define DUMUX_TEST_ANALYTICALSOLVECTORS_HH

namespace Dumux {
/*!
 * \ingroup NavierStokesTests
 * \brief Creates and provides the analytical solution in a form that can be written into the vtk output files
 */
template<class Problem, class Scalar = double>
class NavierStokesAnalyticalSolutionVectors
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using Indices = typename Problem::Indices;

public:
    NavierStokesAnalyticalSolutionVectors(std::shared_ptr<const Problem> problem, Scalar tInitial = 0.0)
    : problem_(problem)
    { update(tInitial); }

    /*!
     * \brief Creates the analytical solution in a form that can be written into the vtk output files
     */
    void update(Scalar time = 0.0)
    {
        analyticalPressure_.resize(problem_->gridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(problem_->gridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(problem_->gridGeometry().numFaceDofs());

        auto fvGeometry = localView(problem_->gridGeometry());

        for (const auto& element : elements(problem_->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                // velocities on faces
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    analyticalVelocityOnFace_[scvf.dofIndex()][scvf.directionIndex()] = problem_->analyticalSolution(scvf.center(), time)[Indices::velocity(scvf.directionIndex())];
                }

                analyticalPressure_[scv.dofIndex()] = problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];

                for (int dirIdx = 0; dirIdx < dimWorld; ++dirIdx)
                    analyticalVelocity_[scv.dofIndex()][dirIdx] = problem_->analyticalSolution(scv.center(), time)[Indices::velocity(dirIdx)];
            }
        }
    }

    /*!
     * \brief Returns the analytical solution for the pressure
     */
    const std::vector<Scalar>& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    const std::vector<VelocityVector>& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    const std::vector<VelocityVector>& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

private:
    std::shared_ptr<const Problem> problem_;

    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;
};
} // end namespace Dumux

#endif
