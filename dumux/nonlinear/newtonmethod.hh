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
#ifndef DUMUX_NEWTONMETHOD_HH
#define DUMUX_NEWTONMETHOD_HH

#include <limits>
#include <dumux/common/exceptions.hh>

#include <dumux/io/vtkmultiwriter.hh>


namespace Dumux
{
/*!
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template <class TypeTag>
class NewtonMethod
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianAssembler)) JacobianAssembler;
public:
    NewtonMethod(Problem &problem)
        : problem_(problem)
    { }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    Problem &problem()
    { return problem_; }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    const Problem &problem() const
    { return problem_; }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    Model &model()
    { return problem().model(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    const Model &model() const
    { return problem().model(); }


    /*!
     * \brief Run the newton method. The controller is responsible
     *        for all the strategic decisions.
     */
    bool execute(NewtonController &ctl)
    {
        try {
            return execute_(ctl);
        }
        catch (const Dune::ISTLError &e) {
            if (ctl.verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            return false;
        }
        catch (const Dumux::NumericalProblem &e) {
            if (ctl.verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            return false;
        };
    }

protected:
    bool execute_(NewtonController &ctl)
    {
        // TODO (?): u shouldn't be hard coded to the model
        SolutionVector &u = model().curSol();
        JacobianAssembler &jacobianAsm = model().jacobianAssembler();

        // tell the controller that we begin solving
        ctl.newtonBegin(*this, u);

        // execute the method as long as the controller thinks
        // that we should do another iteration
        while (ctl.newtonProceed(u))
        {
            // notify the controller that we're about to start
            // a new timestep
            ctl.newtonBeginStep();

            // make the current solution to the old one
            uOld_ = u;

            if (ctl.verbose()) {
                std::cout << "Assembling global jacobian";
                std::cout.flush();
            }
            // linearize the problem at the current solution
            jacobianAsm.assemble(u);

            // solve the resultuing linear equation system
            if (ctl.verbose()) {
                std::cout << "\rSolve Mx = r";
                // Clear the current line using an ansi escape
                // sequence.  for an explaination see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            // set the delta vector to zero before solving the linear system!
            u = 0;

            // ask the controller to solve the linearized system
            ctl.newtonSolveLinear(jacobianAsm.matrix(),
                                  u,
                                  jacobianAsm.residual());

            // update the current solution (i.e. uOld) with the delta
            // (i.e. u). The result is stored in u
            ctl.newtonUpdate(u, uOld_);

            // tell the controller that we're done with this iteration
            ctl.newtonEndStep(u, uOld_);
        }

        // tell the controller that we're done
        ctl.newtonEnd();

        if (!ctl.newtonConverged()) {
            ctl.newtonFail();
            return false;
        }

        ctl.newtonSucceed();
        return true;
    }

private:
    SolutionVector uOld_;
    SolutionVector residual_;

    Problem &problem_;
};

}

#endif
