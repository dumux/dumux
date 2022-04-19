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
#ifndef DUMUX_TEST_FREEFLOW_COMPOSITIONAL_ANALYTICALSOLVECTORS_HH
#define DUMUX_TEST_FREEFLOW_COMPOSITIONAL_ANALYTICALSOLVECTORS_HH

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>

namespace Dumux {
namespace NavierStokesTest {

/*!
 * \ingroup NavierStokesTests
 * \brief Creates and provides the analytical solution in a form that can be written into the vtk output files
 */
template<class MomentumProblem, class MassProblem, class Scalar = double>
class CompositionalAnalyticalSolutionVectors
: public AnalyticalSolutionVectors<MomentumProblem, MassProblem>
{
    using ParentType = AnalyticalSolutionVectors<MomentumProblem, MassProblem>;
    using GridGeometry = std::decay_t<decltype(std::declval<MassProblem>().gridGeometry())>;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
    using MassIndices = typename MassProblem::Indices;
public:
    CompositionalAnalyticalSolutionVectors(std::shared_ptr<const MomentumProblem> momentumProblem,
                                           std::shared_ptr<const MassProblem> massProblem,
                                           Scalar tInitial = 0.0)
    : ParentType(momentumProblem, massProblem, tInitial)
    , massProblem_(massProblem)
    { update(tInitial); }

    /*!
     * \brief Creates the analytical solution in a form that can be written into the vtk output files
     */
    void update(Scalar time = 0.0)
    {
        ParentType::update(time);
        analyticalConcentration_.resize(massProblem_->gridGeometry().numDofs());

        // cell-centers (pressure + velocity)
        {
            auto fvGeometry = localView(massProblem_->gridGeometry());
            const auto gridView = massProblem_->gridGeometry().gridView();
            for (const auto& element : elements(gridView))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    analyticalConcentration_[scv.dofIndex()]
                        = massProblem_->analyticalSolution(scv.dofPosition(), time)[MassIndices::conti0EqIdx + 1];
                }
            }
        }
    }

    /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    const std::vector<Scalar>& analyticalConcentrationSolution() const
    {
        return analyticalConcentration_;
    }

private:
    std::shared_ptr<const MassProblem> massProblem_;

    std::vector<Scalar> analyticalConcentration_;
};

template<class MomentumProblem, class MassProblem>
CompositionalAnalyticalSolutionVectors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>)
-> CompositionalAnalyticalSolutionVectors<MomentumProblem, MassProblem>;

template<class MomentumProblem, class MassProblem, class Scalar>
CompositionalAnalyticalSolutionVectors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>, Scalar)
-> CompositionalAnalyticalSolutionVectors<MomentumProblem, MassProblem, Scalar>;

} // end namespace NavierStokesTest
} // end namespace Dumux

#endif
