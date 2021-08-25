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
template<class Scalar, class FaceSolutionVector, class CellCenterSolutionVector, class GridGeometry, class Problem, class Indices>
class NavierStokesAnalyticalSolutionVectors
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    NavierStokesAnalyticalSolutionVectors(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    { }

    /*!
     * \brief Creates the instationary analytical solution in a form that can be written into the vtk output files
     */
    void update(Scalar time)
    {
        auto instationaryAnalyticalPressureSolution = [&](const SubControlVolume& scv)
        {
            return problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
        };

        auto instationaryAnalyticalVelocitySolution = [&](const SubControlVolumeFace& scvf)
        {
            return problem_->analyticalSolution(scvf.center(), time)[Indices::velocity(scvf.directionIndex())];
        };

        auto instationaryAnalyticalVelocitySolutionAtScvCenter = [&](const SubControlVolume& scv)
        {
            return problem_->analyticalSolution(scv.center(), time);
        };

        createAnalyticalSolution_(instationaryAnalyticalPressureSolution, instationaryAnalyticalVelocitySolution, instationaryAnalyticalVelocitySolutionAtScvCenter);
    }

    /*!
     * \brief Creates the stationary analytical solution in a form that can be written into the vtk output files
     */
    void update()
    {
        auto analyticalPressureSolution = [&](const SubControlVolume& scv)
        {
            return problem_->analyticalSolution(scv.dofPosition())[Indices::pressureIdx];
        };

        auto analyticalVelocitySolution = [&](const SubControlVolumeFace& scvf)
        {
            return problem_->analyticalSolution(scvf.center())[Indices::velocity(scvf.directionIndex())];
        };

        auto analyticalVelocitySolutionAtScvCenter = [&](const SubControlVolume& scv)
        {
            return problem_->analyticalSolution(scv.center());
        };

        createAnalyticalSolution_(analyticalPressureSolution, analyticalVelocitySolution, analyticalVelocitySolutionAtScvCenter);
    }

    /*!
     * \brief Returns the analytical solution for the pressure
     */
    auto& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    auto& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

private:
    /*!
     * \brief Creates the analytical solution in a form that can be written into the vtk output files
     */
    template <class LambdaA, class LambdaB, class LambdaC>
    void createAnalyticalSolution_(const LambdaA& analyticalPressureSolution,
                                  const LambdaB& analyticalVelocitySolution,
                                  const LambdaC& analyticalVelocitySolutionAtScvCenter)
    {
        analyticalPressure_.resize((problem_->gridGeometry()).numCellCenterDofs());
        analyticalVelocity_.resize((problem_->gridGeometry()).numCellCenterDofs());
        analyticalVelocityOnFace_.resize((problem_->gridGeometry()).numFaceDofs());

        for (const auto& element : elements((problem_->gridGeometry()).gridView()))
        {
            auto fvGeometry = localView((problem_->gridGeometry()));
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    analyticalVelocityOnFace_[scvf.dofIndex()][scvf.directionIndex()] = analyticalVelocitySolution(scvf);
                }

                analyticalPressure_[scv.dofIndex()] = analyticalPressureSolution(scv);

                for(int dirIdx = 0; dirIdx < dimWorld; ++dirIdx)
                    analyticalVelocity_[scv.dofIndex()][dirIdx] = analyticalVelocitySolutionAtScvCenter(scv)[Indices::velocity(dirIdx)];
            }
        }
     }

    std::shared_ptr<const Problem> problem_;

    CellCenterSolutionVector analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;
};
} // end namespace Dumux

#endif
