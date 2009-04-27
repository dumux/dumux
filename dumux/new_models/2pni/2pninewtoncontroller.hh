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
class TwoPNINewtonController
    : public NewtonControllerBase<NewtonMethod, TwoPNINewtonController<NewtonMethod> >
{
public:
    typedef TwoPNINewtonController<NewtonMethod>        ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    TwoPNINewtonController(Scalar tolerance = 1e-8,
                               int targetSteps = 9,
                               int maxSteps = 18)
        : ParentType(tolerance, targetSteps, maxSteps)
    {};
};
}

#endif
