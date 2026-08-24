// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \copydoc Dumux::NavierStokesAnalyticalSolutionVectors
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANALYTICALSOLVECTORS_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANALYTICALSOLVECTORS_HH

#include <dune/common/fvector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

namespace Dumux {
/*!
 * \ingroup NavierStokesTests
 * \brief Creates and provides the analytical solution in a form that can be written into the vtk output files
 */
template<class Problem, class Scalar = double>
class NavierStokesAnalyticalSolutionVectors
{
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;
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
        const auto& gridGeometry = Dumux::gridDiscretization(*problem_);

        analyticalPressure_.resize(gridGeometry.numCellCenterDofs());
        analyticalVelocity_.resize(gridGeometry.numCellCenterDofs());
        analyticalVelocityAtDofs_.resize(gridGeometry.numFaceDofs());

        auto fvGeometry = localView(gridGeometry);

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                // velocities on faces
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    analyticalVelocityAtDofs_[scvf.dofIndex()][scvf.directionIndex()] = problem_->analyticalSolution(scvf.center(), time)[Indices::velocity(scvf.directionIndex())];
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
    const std::vector<VelocityVector>& getAnalyticalVelocitySolutionAtDofs() const
    {
        return analyticalVelocityAtDofs_;
    }

private:
    std::shared_ptr<const Problem> problem_;

    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityAtDofs_;
};

template<class Problem>
NavierStokesAnalyticalSolutionVectors(std::shared_ptr<Problem> p)
-> NavierStokesAnalyticalSolutionVectors<Problem>;

template<class Problem, class Scalar>
NavierStokesAnalyticalSolutionVectors(std::shared_ptr<Problem> p, Scalar t)
-> NavierStokesAnalyticalSolutionVectors<Problem, Scalar>;


namespace NavierStokesTest {

/*!
 * \ingroup NavierStokesTests
 * \brief Creates and provides the analytical solution in a form that can be written into the vtk output files
 * \TODO make this the new NavierStokesAnalyticalSolutionVectors once all tests are ported to the new staggered
 */
template<class MomentumProblem, class MassProblem, class Scalar = double>
class AnalyticalSolutionVectors
{
    using MassGridGeometry = typename ProblemTraits<MassProblem>::GridGeometry;
    using MomentumGridGeometry = typename ProblemTraits<MomentumProblem>::GridGeometry;
    static constexpr int dimWorld = MassGridGeometry::GridView::dimensionworld;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using MassIndices = typename MassProblem::Indices;
    using MomIndices = typename MomentumProblem::Indices;
public:
    AnalyticalSolutionVectors(std::shared_ptr<const MomentumProblem> momentumProblem,
                              std::shared_ptr<const MassProblem> massProblem,
                              Scalar tInitial = 0.0)
    : momentumProblem_(momentumProblem)
    , massProblem_(massProblem)
    { update(tInitial); }

    /*!
     * \brief Creates the analytical solution in a form that can be written into the vtk output files
     */
    void update(Scalar time = 0.0)
    {
        const auto& massGridGeometry = Dumux::gridDiscretization(*massProblem_);

        const auto& momentumGridGeometry = Dumux::gridDiscretization(*momentumProblem_);

        analyticalPressure_.resize(massGridGeometry.numDofs());
        analyticalVelocity_.resize(massGridGeometry.gridView().size(0));
        analyticalVelocityAtDofs_.resize(momentumGridGeometry.numDofs());

        // cell-centers (pressure + velocity)
        {
            auto fvGeometry = localView(massGridGeometry);
            const auto gridView = massGridGeometry.gridView();
            for (const auto& element : elements(gridView))
            {
                // output velocity always on elements
                const auto eIdx = massGridGeometry.elementMapper().index(element);
                const auto center = element.geometry().center();
                for (int dirIdx = 0; dirIdx < dimWorld; ++dirIdx)
                        analyticalVelocity_[eIdx][dirIdx]
                            = momentumProblem_->analyticalSolution(center, time)[MomIndices::velocity(dirIdx)];

                // this works for box and cc methods
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                    analyticalPressure_[scv.dofIndex()]
                        = massProblem_->analyticalSolution(scv.dofPosition(), time)[MassIndices::pressureIdx];
            }
        }

        // dof positions (velocity)
        {
            auto fvGeometry = localView(momentumGridGeometry);
            const auto gridView = momentumGridGeometry.gridView();
            for (const auto& element : elements(gridView))
            {
                fvGeometry.bindElement(element);

                if constexpr (MomentumGridGeometry::discMethod == DiscretizationMethods::fcstaggered)
                    for (const auto& scv : scvs(fvGeometry))
                        analyticalVelocityAtDofs_[scv.dofIndex()][scv.dofAxis()]
                            = momentumProblem_->analyticalSolution(scv.center(), time)[MomIndices::velocity(scv.dofAxis())];

                else if constexpr (DiscretizationMethods::isCVFE<typename MomentumGridGeometry::DiscretizationMethod>)
                    for (const auto& scv : scvs(fvGeometry))
                        for (int dirIdx = 0; dirIdx < dimWorld; ++dirIdx)
                            analyticalVelocityAtDofs_[scv.dofIndex()][dirIdx]
                                = momentumProblem_->analyticalSolution(scv.dofPosition(), time)[MomIndices::velocity(dirIdx)];

                else
                    DUNE_THROW(Dune::Exception, "Unknown discretization method: " << MomentumGridGeometry::discMethod);
            }
        }
    }

    /*!
     * \brief Returns the analytical solution for the pressure
     */
    const std::vector<Scalar>& analyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

    /*!
     * \brief Returns the analytical solution for the velocity
     */
    const std::vector<VelocityVector>& analyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

    /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    const std::vector<VelocityVector>& analyticalVelocitySolutionAtDofs() const
    {
        return analyticalVelocityAtDofs_;
    }

private:
    std::shared_ptr<const MomentumProblem> momentumProblem_;
    std::shared_ptr<const MassProblem> massProblem_;

    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityAtDofs_;
};

template<class MomentumProblem, class MassProblem>
AnalyticalSolutionVectors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>)
-> AnalyticalSolutionVectors<MomentumProblem, MassProblem>;

template<class MomentumProblem, class MassProblem, class Scalar>
AnalyticalSolutionVectors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>, Scalar)
-> AnalyticalSolutionVectors<MomentumProblem, MassProblem, Scalar>;


} // end namespace NavierStokesTest
} // end namespace Dumux

#endif
