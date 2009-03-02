//$Id$
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
 * \brief A 2p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2PNI_NEWTON_CONTROLLER_HH
#define DUMUX_2PNI_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/new_newtoncontroller.hh>

namespace Dune {
/*!
 * \brief A 2p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class NewtonMethod>
class TwoPTwoCNINewtonController
    : public NewtonControllerBase<NewtonMethod, TwoPTwoCNINewtonController<NewtonMethod> >
{
public:
    typedef TwoPTwoCNINewtonController<NewtonMethod>            ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;
    typedef typename NewtonMethod::Model                Model;
    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    TwoPTwoCNINewtonController(bool switched = false,
                               Scalar tolerance = 1e-9,
                               int targetSteps = 8,
                               int maxSteps = 12,
                               int switchCount = 0)
        : ParentType(tolerance, targetSteps, maxSteps), switched_(switched), switchCount_(switchCount)
    {};

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // use function of the newtoncontroller
        return ParentType::suggestTimeStepSize(oldTimeStep);
    }

    //! Returns true iff another iteration should be done.
    bool newtonProceed(Function &u)
    {
        // use function of the newtoncontroller
        return ParentType::newtonProceed(u);
    }

    /** \todo Please doc me! */

protected:
    friend class NewtonControllerBase<NewtonMethod, ThisType>;
    //! called by the base class the get an indication of how physical
    //! an iterative solution is 1 means "completely physical", 0 means
    //! "completely unphysical"
    Scalar _physicalness(Function &u)
    {
        return 1.0;
    }
    bool switched_;
    int switchCount_;
};
}

#endif
