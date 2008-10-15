/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
 * \brief Reference implementation of a newton controller solver.
 *
 * Usually for most cases this controller should be sufficient.
 */
#ifndef DUNE_NEWTON_CONTROLLER_HH
#define DUNE_NEWTON_CONTROLLER_HH

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <iostream>
#include <boost/format.hpp>

namespace Dune
{
    /*!
     * \brief Base class for the reference implementation of a newton
     *        controller.
     *
     * If you want to specialise only some methods but are happy with
     * the defaults of the reference controller, derive your
     * controller from this class and simply overload the required
     * methods.
     */
    template <class Model,
              class Implementation>
    class NewtonControllerBase
    {
        typedef typename Model::DomainTraits  DomainTraits;
        typedef typename Model::NewtonTraits  ModelNewtonTraits;

    public:
        typedef typename DomainTraits::Scalar                     Scalar;

        typedef typename ModelNewtonTraits::Function              Function;
        typedef typename ModelNewtonTraits::OperatorAssembler     OperatorAssembler;
        typedef typename OperatorAssembler::RepresentationType    OpAsmRep;

        NewtonControllerBase(Scalar tolerance, // maximum tolerated defect
                             int targetSteps,
                             int maxSteps)
            {
                assert(maxSteps > targetSteps + 3);
                _numSteps = 0;
                _tolerance = tolerance;
                _targetSteps = targetSteps;
                _maxSteps = maxSteps;

                _curPhysicalness = 0;
                _maxPhysicalness = 0;
                _defect = 0;
            };

        //! Returns true iff another iteration should be done.
        bool newtonProceed(Function &u)
            {
                if (_numSteps < 2)
                    return true; // we always do at least two iterations
                else if (_numSteps > _maxSteps)
                    return false; // we have exceeded the allowed number of steps
                else if (newtonConverged())
                    return false; // we have reached the desired defect
                
                // check for the physicalness of the solution
                _curPhysicalness = _asImp()._physicalness(u);
                _curPhysicalness = std::min(_curPhysicalness, 1.0);
                if (_curPhysicalness <= 0)
                    // not physical enough even for a temporary
                    // solution
                    return false;
                else if (_curPhysicalness < ((Scalar) _numSteps)/(_maxSteps - 1)) {
                    // we require that the solution gets more physical
                    // with every step and at the last step the
                    // solution must be completely physical.
                    return false;
                }
                else if (_curPhysicalness < _maxPhysicalness)
                {
                    if (_probationCount > 1) {
                        // an iterative solution was more physical
                        // than the current solution and at least 2
                        // others.
                        return false;
                    }
                    else {
                        // we are physical enough, but some earlier
                        // solution was more physical, so we let the
                        // solver continue on probation.
                        ++_probationCount;
                        return true;
                    }
                }
                else {
                    // everything's fine: the solution is physical
                    // enough for the number of iterations we did and
                    // it is the most physical so far.
                    _maxPhysicalness = _curPhysicalness;
                    _probationCount = std::min(0, _probationCount - 1);

                    return true; // do another round
                };
            }

        //! Returns true iff the defect of the solution is below the
        //! tolerance
        bool newtonConverged()
            {
                return _defect <= _tolerance && _curPhysicalness >= 1.0;
            }

        //! called before the newton method is applied to an equation
        //! system.
        void newtonBegin(Function &u)
            {
                _numSteps = 0;
                _probationCount = 0;
                _maxPhysicalness = 0;
                _oneByMagnitude = 1.0/std::max((*u).two_norm(), 1e-5);
                _defect = 1e100;
            }

        //! indidicates the beginning of a newton iteration
        void newtonBeginStep()
            {
            }

        //! Returns the number of steps done since newtonBegin() was
        //! called
        int newtonNumSteps()
            { return _numSteps; }

        //! Solve the linear equation system Ax - b = 0 for the
        //! current iteration.
        //! Returns true iff the equation system could be solved.
        bool newtonSolveLinear(OperatorAssembler &opAsm,
                               Function &x,
                               Function &b)
            {
                typedef typename Function::RepresentationType           FnRep;
                typedef Dune::MatrixAdapter<OpAsmRep,FnRep,FnRep>       Operator;

                // make make a matrix out of the operator
                Operator A(*opAsm);

                // if the defect of the newton method is large, we do
                // not need to solve the linear approximation
                // accurately. On the other hand, if this is the first
                // newton step, we don't have a meaningful value for the defect
                // yet, so we use the targeted accurracy for the defect.
                Scalar residTol = _tolerance/10;

                // initialize the preconditioner
                Dune::SeqILU0<OpAsmRep,FnRep,FnRep> precond(*opAsm,1.0);
//                Dune::SeqSSOR<OpAsmRep,FnRep,FnRep> precond(*opAsm, 3, 1.0);
//                SeqIdentity<OpAsmRep,FnRep,FnRep> precond(*opAsm);


/*
                for (int i = 0; i < A.getmat().N(); ++i) {
                    for (int k = 0; k < 2; ++k) {
                        for (int j = 0; j < A.getmat().M(); ++j) {
                            for (int l = 0; l < 2; ++l) {
                                if (A.getmat()[i].find(j) != A.getmat()[i].end())
                                    printf("%- #4lf ",
                                           (*A.getmat()[i].find(j))[k][l]);
                                else
                                    printf(".     ");
                            }
                        }
                        printf("\n");
                    }
                }
*/                

                // invert the linear equation system
                Dune::BiCGSTABSolver<FnRep> solver(A, precond, residTol, 10000, 1);
                Dune::InverseOperatorResult result;
                solver.apply(*x, *b, result);
                
                _defect = _oneByMagnitude*((*x).two_norm());
                return result.converged;
            };

        //! Indicates that we're done solving one newton step.
        void newtonEndStep(Function &u, Function &uOld)
            {
                ++_numSteps;
                std::cout << boost::format("Newton iteration %d done: defect=%g, physicalness: %.3f, maxPhysicalness=%.3f\n")
                    %_numSteps%_defect%_curPhysicalness%_maxPhysicalness;
            };

        //! Indicates that we're done solving the equation system.
        void newtonEnd()
            {};

        //! Called when the newton method broke down.
        void newtonFail()
            { 
                _defect = 1e100; _numSteps = _targetSteps*2; 
            }

        //! Suggest a new time stepsize based on the number of newton
        //! iterations required.
        Scalar suggestTimeStepSize(Scalar oldTimeStep) const
            {
                // be agressive reducing the timestep size but
                // conservative when increasing it. the rationale is
                // that we want to avoid failing in the next newton
                // iteration which would require another linerization
                // of the problem.
                if (_numSteps > _targetSteps) {
                    Scalar percent = ((Scalar) _numSteps - _targetSteps)/_targetSteps;
                    return oldTimeStep/(1 + percent);
                }
                else {
                    Scalar percent = ((Scalar) _targetSteps - _numSteps)/_targetSteps;
                    return oldTimeStep*(1 + percent/1.2);
                }
            }
        

    protected:
        // returns the actual implementation for the cotroller we do
        // it this way in order to allow "poor man's virtual methods",
        // i.e. methods of subclasses which can be called by the base
        // class.
        Implementation &_asImp()
            { return *static_cast<Implementation*>(this); }
        const Implementation &_asImp() const
            { return *static_cast<const Implementation*>(this); }

        //! this function is an indication of how "physically
        //! meaningful" a temporary solution is. 0 means it isn't
        //! meaningful at all (maybe because it has highly negative
        //! pressures, etc) and the newton method can be stopped
        //! immediately. Conversly 1 means that the solution is
        //! perfectly physically meaningful (although it doesn't need
        //! to be the solution in any way) and the method can to run.
        //! Values inbetween mean that the funtion is not meaningfull,
        //! but can be tolerated as temporary solution at some
        //! iteration. (The controller assumes that as the method
        //! progresses, the physicallness of the solution must
        //! increase.)
        Scalar _physicalness(Function &u)
            {
                return 1;
            }

        Scalar _defect;
        Scalar _tolerance;

        Scalar _maxPhysicalness;
        Scalar _curPhysicalness;
        Scalar _oneByMagnitude;
        int    _probationCount;
        
        // optimal number of iterations we want to achive
        int    _targetSteps;
        // maximum number of iterations we do before giving up
        int    _maxSteps;
        // actual number of steps done so far
        int    _numSteps;
    };

    //! A reference implementation of a newton method controller
    //!
    //! Basically the only difference to NewtonControllerBase is that
    //! this class can be instanciated more easily.
    template <class Model>
    class NewtonController
        : public NewtonControllerBase<Model, NewtonController<Model> >
    {
    public:
        typedef NewtonController<Model>               ThisType;
        typedef NewtonControllerBase<Model, ThisType> ParentType;

        typedef typename ParentType::Scalar            Scalar;
        typedef typename ParentType::Function          Function;
        typedef typename ParentType::OperatorAssembler OperatorAssembler;

        NewtonController(Scalar tolerance = 1e-5, int maxSteps = 12)
            : ParentType(tolerance, maxSteps)
            {};
    };
}


#endif
