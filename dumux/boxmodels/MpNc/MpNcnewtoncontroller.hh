/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief A MpNc specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_MPNC_NEWTON_CONTROLLER_HH
#define DUMUX_MPNC_NEWTON_CONTROLLER_HH

#include "MpNcproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>
#include <algorithm>

namespace Dumux {

template <class TypeTag, bool enableKinetic /* = false */>
class MpNcNewtonChop
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    enum { numPhases =  GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { fug0Idx = Indices::fug0Idx };
    enum { S0Idx = Indices::S0Idx };
    enum { p0Idx = Indices::p0Idx };

public:
    static void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (int i = 0; i < uLastIter.size(); ++i) {
            for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][S0Idx + phaseIdx],
                                uLastIter[i][S0Idx + phaseIdx]);
            pressureChop_(uCurrentIter[i][p0Idx], uLastIter[i][p0Idx]);
            for (int comp = 0; comp < numComponents; ++comp) {
                pressureChop_(uCurrentIter[i][fug0Idx + comp], uLastIter[i][fug0Idx + comp]);
            }

        }
    };

private:
    static void clampValue_(Scalar &val, Scalar minVal, Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    };

    static void pressureChop_(Scalar &val, Scalar oldVal)
    {
        const Scalar maxDelta = std::max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = std::max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val, Scalar oldVal)
    {
        const Scalar maxDelta = 0.25;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        clampValue_(val, -0.001, 1.001);
    }

};

template <class TypeTag>
class MpNcNewtonChop<TypeTag, /*enableKinetic=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    enum { numPhases =  GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { moleFrac00Idx = Indices::moleFrac00Idx };
    enum { S0Idx = Indices::S0Idx };
    enum { p0Idx = Indices::p0Idx };

public:
    static void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (int i = 0; i < uLastIter.size(); ++i) {
            for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][S0Idx + phaseIdx],
                                uLastIter[i][S0Idx + phaseIdx]);
            pressureChop_(uCurrentIter[i][p0Idx], uLastIter[i][p0Idx]);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    moleFracChop_(uCurrentIter[i][moleFrac00Idx + phaseIdx*numComponents + compIdx],
                                  uLastIter[i][moleFrac00Idx + phaseIdx*numComponents + compIdx]);
                }
            };

        }
    };

private:
    static void clampValue_(Scalar &val, Scalar minVal, Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    };

    static void pressureChop_(Scalar &val, Scalar oldVal)
    {
        const Scalar maxDelta = std::max(oldVal/4.0, 10e3);
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        val = std::max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val, Scalar oldVal)
    {
        const Scalar maxDelta = 0.25;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
        clampValue_(val, -0.001, 1.001);
    }

    static void moleFracChop_(Scalar &val, Scalar oldVal)
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
    typedef MPNCNewtonController<TypeTag> ThisType;
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianAssembler)) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),
        enableEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)),
        enableKinetic = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)),
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)),

        p0Idx = Indices::p0Idx,
        S0Idx = Indices::S0Idx,

        Red = JacobianAssembler::Red
    };

    typedef MpNcNewtonChop<TypeTag, enableKinetic> NewtonChop;

public:
    MPNCNewtonController(const Problem &problem)
        : ParentType(problem)
    { 
        enableChop_ = GET_PARAM(TypeTag, bool, Newton, EnableChop);
        Dune::FMatrixPrecision<>::set_singular_limit(1e-35);
    };

    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        this->writeConvergence_(uLastIter, deltaU);
        this->newtonUpdateRelError(uLastIter, deltaU);

        // compute the vertex and element colors for partial
        // reassembly
        if (enablePartialReassemble) {
            Scalar minReasmTol = 0.5*this->tolerance_;
            Scalar reassembleTol = Dumux::geometricMean(this->error_, minReasmTol);
            reassembleTol = std::max(reassembleTol, minReasmTol);

            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch))) {
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            Scalar lambda = 1.0;
            /*
            if (this->error_ > 10) {
                // Do a "poor man's line search"
                lambda /= std::sqrt((this->numSteps_ + 1)*this->error_/10);
                lambda = std::max(0.2, lambda);
            }
            */

            uCurrentIter = deltaU;
            uCurrentIter *= - lambda;
            uCurrentIter += uLastIter;

            if (this->numSteps_ < 4 && enableChop_) {
                // put crash barriers along the update path at the
                // beginning...
                NewtonChop::chop(uCurrentIter, uLastIter);
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
    };

    bool enableChop_;
};
}

#endif
