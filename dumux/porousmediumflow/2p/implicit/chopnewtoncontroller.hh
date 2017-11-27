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
 * \brief A 2p1cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P_CHOP_NEWTON_CONTROLLER_HH
#define DUMUX_2P_CHOP_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {

/*!
 * \ingroup PNMModel
 * \brief A PNM specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPChopNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
    };

public:
    TwoPChopNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (this->enableShiftCriterion_)
            this->newtonUpdateShift(uLastIter, deltaU);

        this->writeConvergence_(uLastIter, deltaU);

        if (this->useLineSearch_)
        {
            this->lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            for (unsigned int i = 0; i < uLastIter.size(); ++i) {
                uCurrentIter[i] = uLastIter[i];
                uCurrentIter[i] -= deltaU[i];
            }

            if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableChop)) {
                this->chop(uCurrentIter, uLastIter);
            }

            if (this->enableResidualCriterion_)
            {
                SolutionVector tmp(uLastIter);
                this->residual_ = this->method().model().globalResidual(tmp, uCurrentIter);
                this->reduction_ = this->residual_;
                this->reduction_ /= this->initialResidual_;
            }
        }
    }

private:
    void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (unsigned int i = 0; i < uLastIter.size(); ++i) {
            saturationChop_(uCurrentIter[i][saturationIdx],
                            uLastIter[i][saturationIdx]);

            //pressureChop_(uCurrentIter[i][pressureIdx],
            //              uLastIter[i][pressureIdx]);

        }
    };

    void clampValue_(Scalar &val,
                            const Scalar minVal,
                            const Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    };

    void pressureChop_(Scalar &val,
                              const Scalar oldVal)
    {
        const Scalar maxDelta = std::max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = std::max(0.0, val); // don't allow negative pressures
    }

    void saturationChop_(Scalar &val,
                                const Scalar oldVal)
    {
        const Scalar maxDelta = 0.2;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        //clampValue_(val, -0.001, 1.001);
    }

};
}

#endif
