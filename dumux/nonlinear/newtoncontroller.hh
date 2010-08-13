// $Id$
/****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
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
 * \brief Reference implementation of a newton controller solver.
 *
 * Usually for most cases this controller should be sufficient.
 */
#ifndef DUMUX_NEWTON_CONTROLLER_HH
#define DUMUX_NEWTON_CONTROLLER_HH

#include <dumux/common/exceptions.hh>

#include <dune/istl/overlappingschwarz.hh>
//#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include "dune/istl/owneroverlapcopy.hh"

#include <dune/istl/io.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

#include <dumux/common/pardiso.hh>

#include <dumux/common/pdelabpreconditioner.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>


namespace Dumux
{
namespace Properties
{
//! specifies the verbosity of the linear solver (by default it is 0,
//! i.e. it doesn't print anything)
NEW_PROP_TAG(NewtonLinearSolverVerbosity);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

//! specifies whether the update should be done using the line search
//! method instead of the "raw" newton method. whether this property
//! has any effect depends on wether the line search method is
//! implemented for the actual model's newton controller's update()
//! method. By default we do not use line search.
NEW_PROP_TAG(NewtonUseLineSearch);

SET_PROP_DEFAULT(NewtonLinearSolverVerbosity)
{public:
    static const int value = 0;
};

SET_PROP_DEFAULT(NewtonWriteConvergence)
{public:
    static const bool value = false;
};

SET_PROP_DEFAULT(NewtonUseLineSearch)
{public:
    static const bool value = false;
};
};


template <class TypeTag, bool enable>
struct NewtonConvergenceWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    NewtonConvergenceWriter(NewtonController &ctl)
        : ctl_(ctl)
    {
        timeStepIndex_ = 0;
        iteration_ = 0;
        vtkMultiWriter_ = new VtkMultiWriter("convergence");
    }

    ~NewtonConvergenceWriter()
    { delete vtkMultiWriter_; };

    void beginTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    };

    void beginIteration(const GridView &gv)
    {
        ++ iteration_;
        vtkMultiWriter_->beginTimestep(timeStepIndex_ + iteration_ / 100.0,
                                       gv);
    };

    void writeFields(const SolutionVector &uOld,
                     const SolutionVector &deltaU)
    {
        ctl_.method().model().addConvergenceVtkFields(*vtkMultiWriter_, uOld, deltaU);
    };

    void endIteration()
    { vtkMultiWriter_->endTimestep(); };

    void endTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    };

private:
    int timeStepIndex_;
    int iteration_;
    VtkMultiWriter *vtkMultiWriter_;
    NewtonController &ctl_;
};

template <class TypeTag>
struct NewtonConvergenceWriter<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    NewtonConvergenceWriter(NewtonController &ctl)
    {};

    void beginTimestep()
    { };

    void beginIteration(const GridView &gv)
    { };

    void writeFields(const SolutionVector &uOld,
                     const SolutionVector &deltaU)
    { };

    void endIteration()
    { };

    void endTimestep()
    { };
};

/*!
 * \brief The reference implementation of a newton controller.
 *
 * If you want to specialize only some methods but are happy with
 * the defaults of the reference controller, derive your
 * controller from this class and simply overload the required
 * methods.
 */
template <class TypeTag>
class NewtonController
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef NewtonConvergenceWriter<TypeTag, GET_PROP_VALUE(TypeTag, PTAG(NewtonWriteConvergence))>  ConvergenceWriter;

    enum { enableTimeStepRampUp = GET_PROP_VALUE(TypeTag, PTAG(EnableTimeStepRampUp)) };

public:
    NewtonController()
        : endIterMsgStream_(std::ostringstream::out),
          convergenceWriter_(asImp_())
    {
        numSteps_ = 0;
        rampUpSteps_ = 4;

        // maximum acceptable difference of any primary variable
        // between two iterations for convergence
        setRelTolerance(1e-7);
        setTargetSteps(9);
        setMaxSteps(15);
    };

    /*!
     * \brief Set the maximum acceptable difference for convergence of
     *        any primary variable between two iterations
     */
    void setRelTolerance(Scalar tolerance)
    { tolerance_ = tolerance; }

    /*!
     * \brief Set the number of iterations at which the Newton method
     *        should aim at.
     */
    void setTargetSteps(int targetSteps)
    { targetSteps_ = targetSteps; }

    /*!
     * \brief Set the number of iterations after which the Newton
     *        method gives up.
     */
    void setMaxSteps(int maxSteps)
    { maxSteps_ = maxSteps; }
    
    /*!
     * \brief Returns true if another iteration should be done.
     */
    bool newtonProceed(const SolutionVector &u)
    {
        if (numSteps_ < rampUpSteps_)
            return true; // we always do at least two iterations
        else if (numSteps_ >= maxSteps_) {
            // we have exceeded the allowed number of steps.  if the
            // relative error was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of
            // steps
            return error_*4.0 < lastError_ && !asImp_().newtonConverged();
        }
        else if (asImp_().newtonConverged())
            return false; // we are below the desired tolerance

        return true;
    }

    /*!
     * \brief Returns true iff the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
    {
        return error_ <= tolerance_;
    }

    /*!
     * \brief Called before the newton method is applied to an
     *        non-linear system of equations.
     */
    void newtonBegin(NewtonMethod &method, SolutionVector &u)
    {
        method_ = &method;
        numSteps_ = 0;

        if (enableTimeStepRampUp) {
            destTimeStepSize_ = timeManager_().timeStepSize();

            // reduce initial time step size for ramp-up
            timeManager_().setTimeStepSize(destTimeStepSize_ / rampUpSteps_);

            // reduce the tolerance for partial reassembly during the
            // ramp up
            finalTolerance_ = tolerance_;
            tolerance_ = tolerance_*1e8;
        }

        convergenceWriter_.beginTimestep();
    }

    /*!
     * \brief Indidicates the beginning of a newton iteration.
     */
    void newtonBeginStep()
    { lastError_ = error_; }

    /*!
     * \brief Returns the number of steps done since newtonBegin() was
     *        called.
     */
    int newtonNumSteps()
    { return numSteps_; }

    /*!
     * \brief Update the error of the solution compared to the
     *        previous iteration.
     */
    void newtonUpdateRelError(const SolutionVector &uOld,
                              const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        typedef typename SolutionVector::block_type FV;
        error_ = 0;

        int idxI = -1;
        int idxJ = -1;
        for (int i = 0; i < int(uOld.size()); ++i) {
            bool needReassemble = false;
            for (int j = 0; j < FV::size; ++j) {
                Scalar tmp
                    =
                    std::abs(deltaU[i][j])
                    / std::max(std::abs(uOld[i][j]), Scalar(1e-4));
                if (tmp > tolerance_ / 10) {
                    needReassemble = true;
                }
                if (tmp > error_)
                {
                    idxI = i;
                    idxJ = j;
                    error_ = tmp;
                }
            }
            
            model_().jacobianAssembler().markVertexRed(i, 
                                                       needReassemble);
        }

        model_().jacobianAssembler().computeColors();
        error_ = gridView_().comm().max(error_);
    }


    /*!
     * \brief Solve the linear system of equations \f$ \mathbf{A}x - b
     *        = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     */
    template <class Vector>
    void newtonSolveLinear(const JacobianMatrix &A,
                           Vector &u,
                           const Vector &b)
    {
        // if the deflection of the newton method is large, we do not
        // need to solve the linear approximation accurately. Assuming
        // that the initial value for the delta vector u is quite
        // close to the final value, a reduction of 6 orders of
        // magnitude in the defect should be sufficient...
        Scalar residReduction = 1e-6;

        try {
            solveLinear_(A, u, b, residReduction);

            // make sure all processes converged
            int converged = 1;
            gridView_().comm().min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "A process threw MatrixBlockError");
            }
        }
        catch (Dune::MatrixBlockError e) {
            // make sure all processes converged
            int converged = 0;
            gridView_().comm().min(converged);

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
    }

    /*!
     * \brief Update the current solution function with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated here.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behaviour is just to
     * subtract deltaU from uOld.
     *
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     * \param uOld   The solution of the last iteration
     */
    void newtonUpdate(SolutionVector &deltaU, const SolutionVector &uOld)
    {
        writeConvergence_(uOld, deltaU);

        newtonUpdateRelError(uOld, deltaU);
        
        deltaU *= -1;
        deltaU += uOld;
    }


    /*!
     * \brief Indicates that one newton iteration was finished.
     */
    void newtonEndStep(SolutionVector &u, SolutionVector &uOld)
    {
        ++numSteps_;
        
        if (enableTimeStepRampUp) {
            if (numSteps_ < rampUpSteps_) {
                // increase time step size
                Scalar dt = destTimeStepSize_ * (numSteps_ + 1) / rampUpSteps_;
                timeManager_().setTimeStepSize(dt);
            }
            else // if we're not in the ramp-up reduce the target
                 // tolerance
                tolerance_ = finalTolerance_;
        }

        if (verbose())
            std::cout << "\rNewton iteration " << numSteps_ << " done: "
                      << "error=" << error_ << endIterMsg().str() << "\n";
        endIterMsgStream_.str("");
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system of equations.
     */
    void newtonEnd()
    {
        if (enableTimeStepRampUp)
            timeManager_().setTimeStepSize(destTimeStepSize_);
        convergenceWriter_.endTimestep();
    }

    /*!
     * \brief Called if the newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonFail()
    {
        numSteps_ = targetSteps_*2;
    }

    /*!
     * \brief Called when the newton method was sucessful.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed()
    { }

    /*!
     * \brief Suggest a new time stepsize based on the old time step size.
     *
     * The default behaviour is to suggest the old time step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time step.
     */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        Scalar n = numSteps_;
        if (enableTimeStepRampUp) {
            n -= rampUpSteps_/2;
        }
        // be agressive reducing the timestep size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next newton
        // iteration which would require another linearization
        // of the problem.
        if (n > targetSteps_) {
            Scalar percent = (n - targetSteps_)/targetSteps_;
            return oldTimeStep/(1.0 + percent);
        }
        else {
            /*Scalar percent = (Scalar(1))/targetSteps_;
              return oldTimeStep*(1 + percent);
            */
            Scalar percent = (targetSteps_ - n)/targetSteps_;
            return oldTimeStep*(1.0 + percent/1.2);
        }
    }

    /*!
     * \brief Returns a reference to the current newton method
     *        which is controlled by this controller.
     */
    NewtonMethod &method()
    { return *method_; }

    /*!
     * \brief Returns a reference to the current newton method
     *        which is controlled by this controller.
     */
    const NewtonMethod &method() const
    { return *method_; }

    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

    /*!
     * \brief Returns true iff the newton method ought to be chatty.
     */
    bool verbose() const
    { return gridView_().comm().rank() == 0; }

protected:
    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the time manager.
     */
    TimeManager &timeManager_()
    { return problem_().timeManager(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    Model &model_()
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Model &model_() const
    { return problem_().model(); }

    // returns the actual implementation for the cotroller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void writeConvergence_(const SolutionVector &uOld,
                           const SolutionVector &deltaU)
    {
        convergenceWriter_.beginIteration(this->gridView_());
        convergenceWriter_.writeFields(uOld, deltaU);
        convergenceWriter_.endIteration();
    };


    template <class Vector>
    void solveLinear_(const JacobianMatrix &A,
                      Vector &x,
                      const Vector &b,
                      Scalar residReduction)
    {
        int verbosity = GET_PROP_VALUE(TypeTag, PTAG(NewtonLinearSolverVerbosity));
        if (gridView_().comm().rank() != 0)
            verbosity = 0;

#if HAVE_PARDISO
        typedef Dumux::PDELab::ISTLBackend_NoOverlap_Loop_Pardiso<TypeTag> Solver;
        Solver solver(problem_(), 500, verbosity);
#else // !HAVE_PARDISO
#if HAVE_MPI
//        typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GridFunctionSpace> Solver;
//        Solver solver(model_().jacobianAssembler().gridFunctionSpace(), 50000, verbosity);
        typedef Dumux::PDELab::ISTLBackend_NoOverlap_BCGS_ILU<TypeTag> Solver;
        Solver solver(problem_(), 500, verbosity);
#else
        typedef Dumux::PDELab::ISTLBackend_SEQ_BCGS_SSOR Solver;
        Solver solver(500, verbosity);
#endif // HAVE_MPI
#endif // HAVE_PARDISO

        //    Solver solver(model_().jacobianAssembler().gridFunctionSpace(), 500, verbosity);
        Vector bTmp(b);
        solver.apply(A, x, bTmp, residReduction);

        if (!solver.result().converged)
            DUNE_THROW(Dumux::NumericalProblem,
                       "Solving the linear system of equations did not converge.");

        // make sure the solver didn't produce a nan or an inf
        // somewhere. this should never happen but for some strange
        // reason it happens anyway.
        Scalar xNorm2 = x.two_norm2();
        gridView_().comm().sum(xNorm2);
        if (std::isnan(xNorm2) || !std::isfinite(xNorm2))
            DUNE_THROW(Dumux::NumericalProblem,
                       "The linear solver produced a NaN or inf somewhere.");
    }

    std::ostringstream endIterMsgStream_;

    NewtonMethod *method_;

    ConvergenceWriter convergenceWriter_;

    Scalar tolerance_;
    Scalar finalTolerance_;

    Scalar error_;
    Scalar lastError_;

    // number of iterations for the time-step ramp-up
    Scalar rampUpSteps_;
    // time step size after the ramp-up
    Scalar destTimeStepSize_;

    // optimal number of iterations we want to achive
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;
};
} // namespace Dumux

#endif
