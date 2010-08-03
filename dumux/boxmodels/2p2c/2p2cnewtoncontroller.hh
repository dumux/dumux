// $Id: 2p2cnewtoncontroller.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
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
        switchIdx   = Indices::switchIdx
    };

public:
    TwoPTwoCNewtonController()
    {
        this->setRelTolerance(1e-7);
        this->setTargetSteps(9);
        this->setMaxSteps(18);
    };

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    void newtonEndStep(SolutionVector &u, SolutionVector &uOld)
    {
        // call the method of the base class
        this->method().model().updateStaticData(u, uOld);
        ParentType::newtonEndStep(u, uOld);
    }

    void newtonUpdate(SolutionVector &deltaU, const SolutionVector &uOld)
    {
        this->writeConvergence_(uOld, deltaU);
        //Scalar oldRelError = this->error_;
        this->newtonUpdateRelError(uOld, deltaU);

        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch)))
            lineSearchUpdate_(deltaU, uOld);
        else {
            deltaU *= - 1.0;
            deltaU += uOld;
        }
    }


    //! Returns true iff the current solution can be considered to
    //! be acurate enough
    bool newtonConverged()
    {
        if (this->method().model().switched())
            return false;

        return ParentType::newtonConverged();
    };

private:
    void lineSearchUpdate_(SolutionVector &u, const SolutionVector &uOld)
    {
       Scalar lambda = 1.0;
       Scalar globDef;
       SolutionVector tmp(this->method().model(), 0.0);
       Scalar oldGlobDef = this->method().model().globalResidual(tmp);

       int n = 0;
       while (true) {
           u *= -lambda;
           u += uOld;
           globDef = this->method().model().globalResidual(tmp);

           if (globDef < oldGlobDef || lambda <= 1.0/8) {
               this->endIterMsg() << ", defect " << oldGlobDef << "->"  << globDef << "@lambda=2^-" << n;
               return;
           }

           // undo the last iteration
           u -= uOld;
           u /= - lambda;

           // try with a smaller update
           lambda /= 2;
           ++n;
       }
    };
    
};
}

#endif
