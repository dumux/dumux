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
    template <class NewtonMethod, bool lineSearch=false>
    class _NewtonUpdateMethod
    {
        typedef typename NewtonMethod::Model Model;
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
    template <class NewtonMethod>
    class _NewtonUpdateMethod<NewtonMethod, true>
    {
        typedef typename NewtonMethod::Model             Model;
        typedef typename Model::NewtonTraits::Function   Function;
        typedef typename Model::NewtonTraits::Scalar     Scalar;
            
    public:
        _NewtonUpdateMethod(Function &uInitial, Model &model)
            {
//                _oldResidual = ;
            };

            
        bool update(Function &u, 
                    Function &uOld,
                    Model &model)
            {
                // First try a normal newton step
                *u *= - 1;
                *u += *uOld;

                Scalar newGlobalResidual = _computeGlobalResidual(u, model);
                // if the new global residual is larger than the old
                // one, do a line search. this is done by assuming
                // that the square of the residual is a second order
                // polynomial. The at the current newton step the
                // derivative (i.e. the jacobian) and the squared
                // residual are known, and at the next newton
                // iteration, the square of the residual is known.
                // for details, see:
                // J. E. Dennis, R. B. Schnabel: "Numerical methods
                // for unconstrained optimization and nonlinear
                // equations", 1, Prentice-Hall, 1983 pp. 126-127.

                // TODO
                return true;
            };

    private:
        Scalar _oldResidual;
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
        typedef NewtonMethod<Model, useLineSearch>              ThisType;

    public:
        typedef typename JacobianAssembler::RepresentationType  JacobianMatrix;

    public:
        NewtonMethod(Model &model)
            : uOld(model.grid()),
              f(model.grid())
            {
                _residual = NULL;
                _model = NULL;
            }

        ~NewtonMethod()
            {
                delete _residual;
            }

        template <class NewtonController>
        bool execute(Model &model, NewtonController &ctl)
            {
                _model = &model;

                // TODO (?): u shouldn't be hard coded to the model
                Function          &u = model.u();
                LocalJacobian     &localJacobian = model.localJacobian();
                JacobianAssembler &jacobianAsm   = model.jacobianAssembler();
                
                // method to of how updated are done. (either 
                // LineSearch or the plain newton-raphson method)
                Dune::_NewtonUpdateMethod<ThisType, useLineSearch> updateMethod(u, model);;

                _residualUpToDate = false;
                _deflectionTwoNorm = 1e100;
                
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
                        _residualUpToDate = false;
                        _deflectionTwoNorm = (*u).two_norm();
                        
                        // update the current solution. We use either
                        // a line search approach or the plain method.
                        if (!updateMethod.update(u, uOld, model)) {
                            ctl.newtonFail();
                            *model.u() = *model.uOldTimeStep();
                            _model = NULL;
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
                        _model = NULL;
                        return false;
                    }

                    ctl.newtonEndStep(u, uOld);
                }
                // tell the controller that we're done
                ctl.newtonEnd();

                if (!ctl.newtonConverged()) {
                    ctl.newtonFail();
                    _model = NULL;
                    return false;
                }

                _model = NULL;
                return true;
            }

        /*!
         * \brief Returns the current Jacobian matrix.
         */
        const JacobianMatrix &currentJacobian() const
            { return *(_model->jacobianAssembler()); }
        

        /*!
         * \brief Returns the current residual, i.e. the deivation of
         *        the non-linear function from 0 for the current
         *        iteration.
         */
        const Function &residual()
            {
                if (!_residualUpToDate) {
                    if (!_residual)
                        _residual = new Function(f);
                    // update the residual
                    _model->evalGlobalResidual(*_residual);
                    _residualUpToDate = true;
                }
                return *_residual;
            }

        /*!
         * \brief Returns the euclidean norm of the last newton step size.
         */
        Scalar deflectionTwoNorm() const
            { return _deflectionTwoNorm; }

    private:
        Function       uOld;
        Function       f;

        bool          _residualUpToDate;
        Function     *_residual;
        Scalar        _deflectionTwoNorm;
        Model        *_model;
    };
}

#endif
