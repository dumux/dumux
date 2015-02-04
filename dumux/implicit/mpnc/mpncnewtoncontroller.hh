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
 * \brief A MpNc specific controller for the newton solver, which knows
 *       'physically meaningful' solution.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_MPNC_NEWTON_CONTROLLER_HH
#define DUMUX_MPNC_NEWTON_CONTROLLER_HH

#include "mpncproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>
#include <algorithm>

namespace Dumux {

/*!
 * \brief A MpNc specific controller for the newton solver, which knows
 *       'physically meaningful' solution.
 */
template <class TypeTag, bool enableKinetic /* = false */>
class MpNcNewtonChop
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    enum { numPhases =  GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { fug0Idx = Indices::fug0Idx };
    enum { s0Idx = Indices::s0Idx };
    enum { p0Idx = Indices::p0Idx };

public:
    static void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (unsigned int i = 0; i < uLastIter.size(); ++i) {
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][s0Idx + phaseIdx],
                                uLastIter[i][s0Idx + phaseIdx]);
            pressureChop_(uCurrentIter[i][p0Idx], uLastIter[i][p0Idx]);
            for (unsigned int comp = 0; comp < numComponents; ++comp) {
                pressureChop_(uCurrentIter[i][fug0Idx + comp], uLastIter[i][fug0Idx + comp]);
            }

        }
    };

private:
    static void clampValue_(Scalar &val,
                            const Scalar minVal,
                            const Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    };

    static void pressureChop_(Scalar &val,
                              const Scalar oldVal)
    {
        const Scalar maxDelta = std::max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = std::max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val,
                                const Scalar oldVal)
    {
        const Scalar maxDelta = 0.25;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        clampValue_(val, -0.001, 1.001);
    }

};

template <class TypeTag>
class MpNcNewtonChop<TypeTag, /*enableKinetic=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    enum { numPhases =  GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { moleFrac00Idx = Indices::moleFrac00Idx };
    enum { s0Idx = Indices::s0Idx };
    enum { p0Idx = Indices::p0Idx };

public:
    static void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (unsigned int i = 0; i < uLastIter.size(); ++i) {
            for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][s0Idx + phaseIdx],
                                uLastIter[i][s0Idx + phaseIdx]);
            pressureChop_(uCurrentIter[i][p0Idx], uLastIter[i][p0Idx]);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    moleFracChop_(uCurrentIter[i][moleFrac00Idx + phaseIdx*numComponents + compIdx],
                                  uLastIter[i][moleFrac00Idx + phaseIdx*numComponents + compIdx]);
                }
            }

        }
    }

private:
    static void clampValue_(Scalar &val,
                            const Scalar minVal,
                            const Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    };

    static void pressureChop_(Scalar &val,
                              const Scalar oldVal)
    {
        const Scalar maxDelta = std::max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = std::max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val,
                                const Scalar oldVal)
    {
        const Scalar maxDelta = 0.25;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        clampValue_(val, -0.001, 1.001);
    }

    static void moleFracChop_(Scalar &val,
                              const Scalar oldVal)
    {
        // no component mole fraction can change by more than 20% per iteration
        const Scalar maxDelta = 0.20;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        clampValue_(val, -0.001, 1.001);
    }

};

/*!
 * \ingroup Newton
 * \brief A MpNc specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class MPNCNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {    numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {p0Idx = Indices::p0Idx};
    enum {s0Idx = Indices::s0Idx};

    typedef MpNcNewtonChop<TypeTag, enableKinetic> NewtonChop;

public:
    MPNCNewtonController(const Problem &problem)
        : ParentType(problem)
    {
        enableChop_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableChop);
        Dune::FMatrixPrecision<>::set_singular_limit(1e-35);
    };

    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (this->enableShiftCriterion_ || this->enablePartialReassemble_)
            this->newtonUpdateRelError(uLastIter, deltaU);

        // compute the vertex and element colors for partial
        // reassembly
        if (this->enablePartialReassemble_) {
            const Scalar minReasmTol = 1e-2*this->tolerance_;
            const Scalar maxReasmTol = 1e1*this->tolerance_;
            Scalar reassembleTol = std::max(minReasmTol, std::min(maxReasmTol, this->shift_/1e4));

            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        this->writeConvergence_(uLastIter, deltaU);

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, UseLineSearch)) {
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            for (unsigned int i = 0; i < uLastIter.size(); ++i) {
                uCurrentIter[i] = uLastIter[i];
                uCurrentIter[i] -= deltaU[i];
            }

            if (this->numSteps_ < 2 && enableChop_) {
                // put crash barriers along the update path at the
                // beginning...
                NewtonChop::chop(uCurrentIter, uLastIter);
            }

            if (this->enableResidualCriterion_)
            {
                SolutionVector tmp(uLastIter);
                this->reduction_ = this->method().model().globalResidual(tmp, uCurrentIter);
                this->reduction_ /= this->initialResidual_;
            }
        }
    }

private:
    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
       Scalar lambda = 1.0;
       Scalar globDef;

       SolutionVector tmp(uLastIter);
       Scalar oldGlobDef = this->model_().globalResidual(tmp, uLastIter);
       while (true) {
           uCurrentIter  = deltaU;
           uCurrentIter *= -lambda;
           uCurrentIter += uLastIter;

           globDef = this->model_().globalResidual(tmp, uCurrentIter);
           if (globDef < oldGlobDef || lambda <= 1.0/64) {
               this->endIterMsg() << ", defect " << oldGlobDef << "->"  << globDef << "@lambda=1/" << 1/lambda;
               return;
           }

           // try with a smaller update
           lambda /= 2;
       }
    }

    bool enableChop_;
};
}

#endif
