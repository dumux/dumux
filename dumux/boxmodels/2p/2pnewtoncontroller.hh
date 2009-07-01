// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief A newton controller for two-phase problems.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P_NEWTON_CONTROLLER_HH
#define DUMUX_2P_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dune {
/*!
 * \ingroup TwoPBoxModel
 * \brief A newton controller for two-phase problems.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class NewtonMethod>
class TwoPNewtonController
    : public NewtonControllerBase<NewtonMethod, TwoPNewtonController<NewtonMethod> >
{
public:
    typedef TwoPNewtonController<NewtonMethod>            ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    TwoPNewtonController(Scalar tolerance = 1e-5,
                         int targetSteps = 8,
                         int maxSteps = 12)
        : ParentType(tolerance, targetSteps, maxSteps)
    {};

protected:
    friend class NewtonControllerBase<NewtonMethod, ThisType>;
    //! called by the base class the get an indication of how physical
    //! an iterative solution is 1 means "completely physical", 0 means
    //! "completely unphysical"
    Scalar physicalness_(Function &u)
    {
        const Scalar satNormFactor = 2.5;

        // the maximum distance of a saturation value to a physically
        // meaningful value.
        Scalar maxSatDelta = 0;
        Scalar sat;
        //                Scalar pressure;

        for (int idx = 0; idx < (int) (*u).size(); idx++)  {
            // TODO: Don't expect the saturation at the second
            // position in the primary var vector
            // pressure = (*u)[idx][0];
            sat = (*u)[idx][1];

            if (sat < 0) {
                maxSatDelta = std::max(maxSatDelta, std::abs(sat));
            }
            else if (sat > 1) {
                maxSatDelta = std::max(maxSatDelta, std::abs(sat - 1));
            }

            // (so far we ignore the phase pressure)
        }

        // we accept solutions up to 0.2 percent bigger than 1
        // or smaller than 0 as being physical for numerical
        // reasons...
        Scalar phys = 1.002 - maxSatDelta/satNormFactor;

        // we never return exactly zero, since we want to
        // allow solutions which are "very close" to a
        // physically meaningful one
        return std::min(1.0, phys);
    }
};
}

#endif
