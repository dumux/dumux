// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
 * \brief A 2p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P2C_NEWTON_CONTROLLER_HH
#define DUMUX_2P2C_NEWTON_CONTROLLER_HH

#include "2p2cproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup TwoPTwoCModel
 * \brief A 2p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPTwoCNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    enum { enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)) };

public:
    TwoPTwoCNewtonController()
    {
        this->setRelTolerance(1e-7);
        this->setTargetSteps(9);
        this->setMaxSteps(18);
    };


    /*!
     * \brief
     * Suggest a new time step size based either on the number of newton
     * iterations required or on the variable switch
     *
     * \param u The current global solution vector
     * \param uLastIter The previous global solution vector
     *
     */
    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        // call the method of the base class
        this->method().model().updateStaticData(uCurrentIter, uLastIter);
        ParentType::newtonEndStep(uCurrentIter, uLastIter);
    }

    /*!
     * \brief Update the current solution function with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated here.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behaviour is just to
     * subtract deltaU from uLastIter.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution after the last Newton iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        this->writeConvergence_(uLastIter, deltaU);

        this->newtonUpdateRelError(uLastIter, deltaU);

        // compute the vertex and element colors for partial
        // reassembly
        if (enablePartialReassemble) {
            Scalar reassembleTol = Dumux::geometricMean(this->error_, 0.1*this->tolerance_);
            reassembleTol = std::max(reassembleTol, 0.1*this->tolerance_);
            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch)))
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        else {
            uCurrentIter = uLastIter;
            uCurrentIter -= deltaU;
        }
    }


    /*!
     * \brief
     * Returns true if the current solution can be considered to
     * be accurate enough
     */
    bool newtonConverged()
    {
        if (this->method().model().switched())
            return false;

        return ParentType::newtonConverged();
    };

private:
    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
       Scalar lambda = 1.0;
       Scalar globDef;

       // calculate the residual of the current solution
       SolutionVector tmp(uLastIter);
       Scalar oldGlobDef = this->method().model().globalResidual(tmp, uLastIter);

       while (true) {
           uCurrentIter = deltaU;
           uCurrentIter *= -lambda;
           uCurrentIter += uLastIter;

           // calculate the residual of the current solution
           globDef = this->method().model().globalResidual(tmp, uCurrentIter);

           if (globDef < oldGlobDef || lambda <= 1.0/8) {
               this->endIterMsg() << ", defect " << oldGlobDef << "->"  << globDef << "@lambda=" << lambda;
               return;
           }

           // try with a smaller update
           lambda /= 2;
       }
    };

};
}

#endif
