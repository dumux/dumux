// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 *
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

#ifdef DUNE_NEWTONMETHOD_DEPRECATED_HH
# error "Never use the old and the new newton method in one program!"
#endif

#include <limits>
#include <dumux/exceptions.hh>


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
    _NewtonUpdateMethod(NewtonMethod &newton,
                        Function &uInitial,
                        Model &model)
    {};

    bool update(NewtonMethod &newton,
                Function &u,
                Function &uOld,
                Model &model)
    {
        *u *= -1.0;
        *u += *uOld;

        // invalidate the current residual vector
        newton.setResidualObsolete();
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
public:
    typedef typename NewtonMethod::Model             Model;
    typedef typename Model::NewtonTraits::Function   Function;
    typedef typename Model::NewtonTraits::Scalar     Scalar;

    _NewtonUpdateMethod(NewtonMethod &newton,
                        Function &uInitial,
                        Model &model)
    {
        uNormMin_ = std::numeric_limits<Scalar>::max();
    };


    bool update(NewtonMethod &newton,
                Function &u,
                Function &uOld,
                Model &model)
    {
        Scalar lambda = 1.0;

        Scalar uNorm = (*u).two_norm();
        if (uNorm < uNormMin_) {
            uNormMin_ = uNorm;
            
            lambda = 1.0;
        }
        else { 
            lambda = uNormMin_ / uNorm;
            if (lambda < 0.05) {
                DUNE_THROW(NumericalProblem, 
                           "Increase of relative defect of " << 1./lambda << " is too much!");
            };
            std::cout << boost::format("Newton: use line search, lambda=%f\n")%lambda;
        }
        
        *u *= - lambda;
        *u += *uOld;

        return true;
    };

private:
    Scalar uNormMin_;
};

/*!
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template<class ModelT, bool useLineSearch=true>
class NewtonMethod
{
public:
    typedef ModelT Model;

    //    private:
    typedef typename Model::NewtonTraits::Function          Function;
    typedef typename Model::NewtonTraits::LocalJacobian     LocalJacobian;
    typedef typename Model::NewtonTraits::JacobianAssembler JacobianAssembler;
    typedef typename Model::NewtonTraits::Scalar            Scalar;
    typedef NewtonMethod<Model, useLineSearch>           ThisType;

public:
    typedef typename JacobianAssembler::RepresentationType  JacobianMatrix;

public:
    NewtonMethod(Model &model)
        : uOld(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0),
          f(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0)
    {
        deflectionTwoNorm_ = 1e100;
        residual_ = NULL;
        model_ = NULL;
    }

    ~NewtonMethod()
    {
        delete residual_;
    }

    /*!
     * \brief Returns a reference to the current numeric model.
     */
    Model &model()
    { return *model_; }

    /*!
     * \brief Returns a reference to the current numeric model.
     */
    const Model &model() const
    { return *model_; }

    /*!
     * \brief Returns true iff the newton method should be verbose
     *        about what it is going on;
     */
    bool verbose() const
    { return (model().gridView().comm().rank() == 0); }

    /*!
     * \brief Run the newton method. The controller is responsible
     *        for all the strategic decisions.
     */
    template <class NewtonController>
    bool execute(Model &model, NewtonController &ctl)
    {
        try {
            return execute_(model, ctl);
        }
        catch (const Dune::NumericalProblem &e) {
            if (verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            model_ = NULL;
            return false;
        };
    };

    /*!
     * \brief Returns the current Jacobian matrix.
     */
    const JacobianMatrix &currentJacobian() const
    { return *(model_->jacobianAssembler()); }


    /*!
     * \brief This method causes the residual to be recalcuated
     *        next time the residual() method is called. It is
     *        internal and only used by some update methods such
     *        as LineSearch.
     */
    void setResidualObsolete(bool yesno=true)
    { residualUpToDate_ = !yesno; };

    /*!
     * \brief Returns the current residual, i.e. the deivation of
     *        the non-linear function from 0 for the current
     *        iteration.
     */
    Function &residual()
    {
        if (!residualUpToDate_) {
            if (!residual_)
                residual_ = new Function(model_->grid(), 
                                         model().grid().overlapSize(0) == 0);
            // update the residual
            model_->evalGlobalResidual(*residual_);
            residualUpToDate_ = true;
        }

        return *residual_;
    }

    /*!
     * \brief Returns the euclidean norm of the last newton step size.
     */
    Scalar deflectionTwoNorm() const
    { return deflectionTwoNorm_; }

protected:
    template <class NewtonController>
    bool execute_(Model &model, NewtonController &ctl)
    {
        model_ = &model;

        // TODO (?): u shouldn't be hard coded to the model
        Function          &u             = model.curSolFunction();
        LocalJacobian     &localJacobian = model.localJacobian();
        JacobianAssembler &jacobianAsm   = model.jacobianAssembler();

        // method to of how updated are done. (either
        // LineSearch or the plain newton-raphson method)
        Dune::_NewtonUpdateMethod<ThisType, useLineSearch> updateMethod(*this, u, model);;

        residualUpToDate_ = false;
        deflectionTwoNorm_ = 1e100;

        // tell the controller that we begin solving
        ctl.newtonBegin(this, u);

        // execute the method as long as the controller thinks
        // that we should do another iteration
        while (ctl.newtonProceed(u))
        {
            // notify the controller that we're about to start
            // a new timestep
            ctl.newtonBeginStep();

            // make the current solution to the old one
            *uOld = *u;
            *f = 0;

            // linearize the problem at the current solution
            jacobianAsm.assemble(localJacobian, u, f);

            // solve the resultuing linear equation system
            ctl.newtonSolveLinear(*jacobianAsm, u, *f);

            Scalar tmp = (*u).two_norm2();
            tmp = model.gridView().comm().sum(tmp);
            deflectionTwoNorm_ = sqrt(tmp);

            // update the current solution. We use either
            // a line search approach or the plain method.
            updateMethod.update(*this, u, uOld, model);

            ctl.newtonEndStep(u, uOld);
        }
        // tell the controller that we're done
        ctl.newtonEnd();

        if (!ctl.newtonConverged()) {
            ctl.newtonFail();
            model_ = NULL;
            return false;
        }

        model_ = NULL;
        return true;
    }


private:
    Function       uOld;
    Function       f;

    bool          residualUpToDate_;
    Function     *residual_;
    Scalar        deflectionTwoNorm_;
    Model        *model_;
    bool          verbose_;
};
}

#endif
