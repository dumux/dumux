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
#ifndef DUMUX_2P2C_NEWTON_CONTROLLER_HH
#define DUMUX_2P2C_NEWTON_CONTROLLER_HH

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
class TwoPTwoCNewtonController
    : public NewtonControllerBase<NewtonMethod, TwoPTwoCNewtonController<NewtonMethod> >
{
    typedef typename NewtonMethod::Model Model;
    typedef typename Model::TwoPTwoCTraits TwoPTwoCTraits;

    enum {
        pressureIdx = TwoPTwoCTraits::pressureIdx,
        switchIdx   = TwoPTwoCTraits::switchIdx
    };

public:
    typedef TwoPTwoCNewtonController<NewtonMethod>        ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    TwoPTwoCNewtonController(Scalar tolerance = 1e-7,
                             int targetSteps = 9,
                             int maxSteps = 18)
        : ParentType(tolerance, targetSteps, maxSteps)
    {};

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    void newtonEndStep(Function &u, Function &uOld)
    {
        // call the method of the base class
        ParentType::model().localJacobian().updateStaticData(u, uOld);
        ParentType::newtonEndStep(u, uOld);
    }

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // use function of the newtoncontroller
        return ParentType::suggestTimeStepSize(oldTimeStep);
    }

    //! Returns true iff the current solution can be considered to
    //! be acurate enough
    bool newtonConverged()
    {
        if (ParentType::model().switched())
            return false;

        return ParentType::newtonConverged();
    };

    //! Returns true iff another iteration should be done.
    bool newtonProceed(Function &u)
    {
        return ParentType::newtonProceed(u);
    }

    /** \todo Please doc me! */

protected:
    friend class NewtonControllerBase<NewtonMethod, ThisType>;
    //! called by the base class the get an indication of how physical
    //! an iterative solution is 1 means "completely physical", 0 means
    //! "completely unphysical"
    Scalar physicalness_(Function &u)
    {
        // the maximum distance of a Sn value to a physically
        // meaningful value.
        Scalar maxSwitchVarDelta = 0;
        Scalar maxPwDelta = 0;

        for (int idx = 0; idx < (int) (*u).size(); idx++)
        {
            Scalar pressure = (*u)[idx][pressureIdx];
            Scalar switchVar = (*u)[idx][switchIdx];

            if (switchVar < 0.0) {
                maxSwitchVarDelta = std::max(maxSwitchVarDelta,
                                             std::abs(switchVar));
            }
            else if (switchVar > 1.0) {
                maxSwitchVarDelta = std::max(maxSwitchVarDelta,
                                             std::abs(switchVar - 1));
            }

            if (pressure < 0.0){
                maxPwDelta = std::max(maxPwDelta,
                                      std::abs(pressure));
            }
        }

        // we accept solutions up to 1 percent bigger than 1
        // or smaller than 0 as being physical for numerical
        // reasons...
        Scalar phys = 1.0001 - maxSwitchVarDelta*100 - maxPwDelta/1e5;

        return std::min(1.0, phys);
    }
};
}

#endif
