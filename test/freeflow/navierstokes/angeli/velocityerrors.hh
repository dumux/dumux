// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Errors of the velocity gradient and divergence for the face-centered staggered scheme on a Cartesian grid.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANGELI_VELOCITY_ERRORS_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANGELI_VELOCITY_ERRORS_HH

#include <cmath>
#include <fstream>
#include <memory>
#include <string>

#include <dumux/io/format.hh>
#include <dumux/common/typetraits/problem.hh>

namespace Dumux::NavierStokesTest {

/*!
 * \ingroup NavierStokesTests
 * \brief Discrete gradient and divergence errors of a staggered velocity field.
 *
 * The gradient is evaluated with the difference quotients of the viscous flux: \f$ \partial u_i / \partial x_i \f$
 * in the cells between opposite faces, \f$ \partial u_i / \partial x_j \f$ between neighboring parallel faces or a
 * face and the boundary, where the analytical solution is used. Each quotient is compared with the
 * analytical gradient at the midpoint of the two positions and weighted with the area it covers.
 * The problem has to provide `analyticalVelocityGradient(globalPos, time)`.
 */
template<class MomentumProblem>
class VelocityGradientErrors
{
    using GridGeometry = typename ProblemTraits<MomentumProblem>::GridGeometry;
    using Scalar = typename GridGeometry::GridView::ctype;

public:
    VelocityGradientErrors(std::shared_ptr<const MomentumProblem> problem, const std::string& name)
    : problem_(problem)
    , fileName_(name + "_velocity_error.csv")
    {
        std::ofstream file(fileName_, std::ios::trunc);
        file << "time, gradientError2, gradientReference2, divergence2\n";
    }

    //! Compute the squared norms at the given time and append them to the csv file
    template<class MomentumSolution>
    void update(const MomentumSolution& x, Scalar time)
    {
        Scalar gradientError2 = 0.0;
        Scalar gradientReference2 = 0.0;
        Scalar divergence2 = 0.0;

        const auto& gridGeometry = problem_->gridGeometry();
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            Scalar divergence = 0.0;

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto axis = scv.dofAxis();
                for (const auto& scvf : scvfs(fvGeometry, scv))
                {
                    if (scvf.isFrontal() && scvf.boundary())
                        continue;

                    const auto normalAxis = scvf.normalAxis();
                    auto outsidePos = scv.dofPosition();
                    Scalar outsideValue = 0.0;
                    if (scvf.boundary())
                    {
                        outsidePos[normalAxis] = scvf.center()[normalAxis];
                        outsideValue = problem_->analyticalSolution(outsidePos, time)[axis];
                    }
                    else
                    {
                        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                        // the two halves of the control volume of a face belong to the same dof
                        if (outsideScv.dofIndex() == scv.dofIndex())
                            continue;

                        outsidePos = outsideScv.dofPosition();
                        outsideValue = x[outsideScv.dofIndex()][0];
                    }

                    const Scalar distance = outsidePos[normalAxis] - scv.dofPosition()[normalAxis];
                    const Scalar gradient = (outsideValue - x[scv.dofIndex()][0])/distance;

                    auto midPoint = outsidePos;
                    midPoint += scv.dofPosition();
                    midPoint *= 0.5;
                    const Scalar exact = problem_->analyticalVelocityGradient(midPoint, time)[axis][normalAxis];

                    // inner quotients are visited from both sides
                    const Scalar weight = (scvf.boundary() ? 1.0 : 0.5)*scvf.area()*std::abs(distance);
                    gradientError2 += weight*(gradient - exact)*(gradient - exact);
                    gradientReference2 += weight*exact*exact;

                    if (scvf.isFrontal())
                        divergence += 0.5*gradient;
                }
            }

            divergence2 += element.geometry().volume()*divergence*divergence;
        }

        std::ofstream file(fileName_, std::ios::app);
        file << Fmt::format("{:.5e}, {:.12e}, {:.12e}, {:.12e}\n", time, gradientError2, gradientReference2, divergence2);
    }

private:
    std::shared_ptr<const MomentumProblem> problem_;
    std::string fileName_;
};

} // end namespace Dumux::NavierStokesTest

#endif
