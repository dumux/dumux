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
#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

namespace Dune
{
    /*!
     * \internal
     * \brief This internal class updates the current iteration u_k of
     *        the solution of the newton method to u_{k+1}. It is not
     *        in the actual execute method because we need
     *        specialization in order not to call interfaces which
     *        might not be present in the model if we do not use the
     *        line search newton method.
     *
     * This is the version without line search.
     */
    template <class Model, bool lineSearch=false>
    class _NewtonUpdateMethod
    {
        typedef typename Model::NewtonTraits::Function   Function;
        typedef typename Model::NewtonTraits::Scalar     Scalar;

    public:
        _NewtonUpdateMethod(Function &uInitial, Model &model)
            {};
            
        bool update(Function &u, Function &uOld, Model &model)
            {
                *u *= -1.0;
                *u += *uOld;
                return true;
            };
    };

    /*!
     * \internal
     * \brief This internal class updates the current iteration u_k of
     *        the solution of the newton method to u_{k+1}. It is not
     *        in the actual execute method because we need
     *        specialization in order not to call interfaces which
     *        might not be present in the model if we do not use the
     *        line search newton method.
     *
     * This is the version with line search. Models using this require
     * an evalGlobalDefect() method.
     */
    template <class Model>
    class _NewtonUpdateMethod<Model, true>
    {
        typedef typename Model::NewtonTraits::Function   Function;
        typedef typename Model::NewtonTraits::Scalar     Scalar;
            
    public:
        _NewtonUpdateMethod(Function &uInitial, Model &model)
            : _globalDefect(model.grid())
            {
                _lambda = 1.0;
                _curGlobalResidual = _computeGlobalResidual(uInitial, model);
            };

            
        bool update(Function &u, 
                    Function &uOld,
                    Model &model)
            {                
/*
                if (recursionDepth == 0) {
                    // if this is the outermost call to this function
                    // first try to chose step size 1.0 which
                    // corrosponds to the plain newton-raphson method
                    _lambda = 1.0;
                }
                else if (recursionDepth > 3) {
                    // if the step size gets too small, cancel the
                    // update
                    return true;
                }
*/

                // do the actual update
                *u *= -_lambda;
                *u += *uOld;

                if (_lambda < 1.0/32) {
                    // if the step size gets too small, cancel the
                    // update
                    return true;
                }

                Scalar newGlobalResidual = _computeGlobalResidual(u, model);
                if (newGlobalResidual > _curGlobalResidual) {
                    // if the square norm of the new global residual
                    // is bigger than the one of the last iteration,
                    // use a smaller step size lambda for the next
                    // iteration.
                    *u = *uOld;
                    _lambda /= 2;
                }
                
                _curGlobalResidual = newGlobalResidual;
                std::cout << boost::format("Newton line search lambda: %f\n")%_lambda;
                return true;
            };

    protected:
        Scalar _computeGlobalResidual(Function &u, Model &model)
            {
                // calculate global residual of a solution
                model.evalGlobalDefect(_globalDefect);
                Scalar globalResidual = (*_globalDefect).two_norm()/2.0;
//                globalResidual *= globalResidual;
                model.grid().comm().sum(&globalResidual, 1);
                std::cout << boost::format("Newton: globalResidual: %f\n")%globalResidual;
//                residualWeight = 1.0/std::max(globalResidual, (Scalar) 1.0e-8);
//                globalResidual *= residualWeight;
                
                return globalResidual;                
            }

        Scalar    _curGlobalResidual;
        Scalar    _lambda; // scaling factor of the step for the update
        Function  _globalDefect; // temporary variable for _computeGlobalResidual
    };

    /*!
     * \brief The algorithmic part of the multi dimensional newton method.
     *
     * In order to use the method you need a \ref NewtonController.
     */
    template<class Model, class NewtonController, bool useLineSearch=true>
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
                
                // method to of how updated are done. (either 
                // LineSearch or the plain newton-raphson method)
                Dune::_NewtonUpdateMethod<Model, useLineSearch> updateMethod(u, model);;

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
                        // update the current solution. We use either
                        // a line search approach or the plain method.
                        if (!updateMethod.update(u, uOld, model)) {
                            ctl.newtonFail();
                            *model.u() = *model.uOldTimeStep();
                            return false;
                        }
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
