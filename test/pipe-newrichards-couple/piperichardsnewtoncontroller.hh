/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan, Andreas Lauser                        *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: onur.dogan _at_ iws.uni-stuttgart.de                             *
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
 * \brief A newton controller for the Richards model which is coupled
 *        with a pipe flow.
 */
#ifndef DUMUX_PIPE_RICHARDS_NEWTON_CONTROLLER_HH
#define DUMUX_PIPE_RICHARDS_NEWTON_CONTROLLER_HH

#include <algorithm>
#include <dumux/nonlinear/new_newtoncontroller.hh>

namespace Dune {
/*!
 * \brief A newton controller for the Richards model which is coupled
 *        with a pipe flow.
 */
template <class NewtonMethod>
class PipeRichardsNewtonController
    : public NewtonControllerBase<NewtonMethod,
                                  PipeRichardsNewtonController<NewtonMethod> >
{
    typedef PipeRichardsNewtonController<NewtonMethod>     ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>   ParentType;

public:
    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    PipeRichardsNewtonController(Scalar tolerance = 1e-5,
                                 int targetSteps = 8,
                                 int maxSteps = 12,
                                 Scalar maxStepSize=1e100)
        : ParentType(tolerance, targetSteps, maxSteps),
          maxStepSize_(maxStepSize)
    {};


    //! Suggest a new time stepsize based on the number of newton
    //! iterations required for the last time step and the old time
    //! step size.
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        Scalar tmp = ParentType::suggestTimeStepSize(oldTimeStep);
        return std::min(tmp, maxStepSize_);
    }


    //! Indicates that we're done solving one newton step.
    void newtonBeginStep()
    {
        // update the pipe flow before each newton-raphson
        // iteration
        this->model().problem().updatePipeFlow();

        // recycle the base class' beginStep method
        ParentType::newtonBeginStep();
    };
private:
    Scalar maxStepSize_;
};
}

#endif
