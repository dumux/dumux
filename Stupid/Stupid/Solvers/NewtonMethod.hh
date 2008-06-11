/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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
 *
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a \ref NewtonController.
 */
#ifndef STUPID_NEWTONMETHOD_HH
#define STUPID_NEWTONMETHOD_HH

namespace Stupid
{

    /*!
     * \brief The algorithmic part of the multi dimensional newton method.
     *
     * In order to use the method you need a \ref NewtonController.
     */
    template<class Model, class NewtonController>
    class NewtonMethod
    {
        typedef typename Model::NewtonTraits::Function          Function;
        typedef typename Model::NewtonTraits::OperatorAssembler OperatorAssembler;
        typedef typename Model::NewtonTraits::LocalJacobian     LocalJacobian;


    public:
        NewtonMethod(Model &model)
            : uOld(model.grid()),
              f(model.grid())
            {
            }

    private:
        Function       uOld;
        Function       f;

    public:
        bool execute(Model &model, NewtonController &ctl)
            {
                // TODO (?): u shouldn't be hard coded to the model
                Function      &u = model.u();
                LocalJacobian &localJacobian = model.localJacobian();
                OperatorAssembler &opAsm = model.opAsm();

                // tell the controller that we begin solving
                ctl.newtonBegin(u);

                // execute the method as long as the controller thinks
                // that we should do another iteration
                while (ctl.newtonProceed(u))
                {
                    // notify the controller that we're about to start
                    // a new timestep
                    ctl.newtonBeginStep();

                    // make the current soltion to the old one
                    *uOld = *u;
                    *f = 0;

                    // linearize the problem at the current solution
                    localJacobian.clearVisited();
                    opAsm.assemble(localJacobian, u, f);

                    // solve the resultuing linear equation system
                    if (ctl.newtonSolveLinear(opAsm, u, f)) {
                        // update the current solution
                        *u *= -1.0;
                        *u += *uOld;
                    }
                    else {
                        // couldn't solve the current linearization
                        ctl.newtonFail();
                        // reset the current solution of the model to the
                        // solution of the last time step in order not to
                        // spoil a possible rerun of the newton method
                        // by artifacts of the current run. (the solution
                        // which didn't converge might be completely unphysical,
                        // but we would like to start from something which is
                        // physical.)
                        *model.u() = *model.uOldTimeStep();
                        return false;
                    }

                    ctl.newtonEndStep(u, uOld);
                }
                // tell the controller that we're done
                ctl.newtonEnd();

                if (!ctl.newtonConverged()) {
                    ctl.newtonFail();
                    return false;
                }
                
                return true;
            }
    };
}

#endif
