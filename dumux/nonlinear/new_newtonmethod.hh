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

                if (_lambda < 1.0/8) {
                    // if the step size gets too small, cancel the
                    // update
                    return true;
                }

                Scalar newGlobalResidual = _computeGlobalResidual(u, model);
                std::cout << boost::format("Newton lambda: %f, residual ratios: %f\n")%_lambda%(newGlobalResidual/_curGlobalResidual);
                if (newGlobalResidual > _curGlobalResidual*1.05) {
                    // if the square norm of the new global residual
                    // is bigger than the one of the last iteration,
                    // use a smaller step size lambda for the next
                    // iteration.
                    _lambda /= 2;
                }
                else if (newGlobalResidual < _curGlobalResidual/1.5) {
                    _lambda = std::min(Scalar(1.0), _lambda*2);
                }
                
                _curGlobalResidual = newGlobalResidual;
                return true;
            };

    protected:
        Scalar _computeGlobalResidual(Function &u, Model &model)
            {
                // calculate global residual of a solution
                model.evalGlobalDefect(_globalDefect);
                Scalar globalResidual = (*_globalDefect).two_norm2();
//                globalResidual *= globalResidual;
                model.grid().comm().sum(&globalResidual, 1);
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
    template<class ModelT, bool useLineSearch=false>
    class NewtonMethod
    {
    public:
        typedef ModelT Model;

    private:
        typedef typename Model::NewtonTraits::Function          Function;
        typedef typename Model::NewtonTraits::LocalJacobian     LocalJacobian;
        typedef typename Model::NewtonTraits::JacobianAssembler JacobianAssembler;
        typedef typename Model::NewtonTraits::Scalar            Scalar;

    public:
        NewtonMethod(Model &model)
            : uOld(model.grid()),
              f(model.grid())
            {
            }

        template <class NewtonController>
        bool execute(Model &model, NewtonController &ctl)
            {
                _defect = 1e100;

                // TODO (?): u shouldn't be hard coded to the model
                Function          &u = model.u();
                LocalJacobian     &localJacobian = model.localJacobian();
                JacobianAssembler &jacobianAsm   = model.jacobianAssembler();
                
                // method to of how updated are done. (either 
                // LineSearch or the plain newton-raphson method)
                Dune::_NewtonUpdateMethod<Model, useLineSearch> updateMethod(u, model);;

                // tell the controller that we begin solving
                ctl.newtonBegin(this, u);

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
                    jacobianAsm.assemble(localJacobian, u, f);

                    // solve the resultuing linear equation system
                    if (ctl.newtonSolveLinear(*jacobianAsm, *u, *f)) {
                        _defect = (*u).two_norm();

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
        
        Scalar defect() const
            { return _defect; }

    private:
        Function       uOld;
        Function       f;
        Scalar        _defect;
    };
}

#endif
