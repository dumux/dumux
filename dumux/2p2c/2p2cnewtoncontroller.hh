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
    public:
        typedef TwoPTwoCNewtonController<NewtonMethod>        ThisType;
        typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

        typedef typename NewtonMethod::Model           Model;
        typedef typename ParentType::Scalar            Scalar;
        typedef typename ParentType::Function          Function;
        typedef typename ParentType::JacobianAssembler JacobianAssembler;

        TwoPTwoCNewtonController(Scalar tolerance = 1e-5,
                                 int targetSteps = 8,
                                 int maxSteps = 12)
            : ParentType(tolerance, targetSteps, maxSteps), minStepsSwitch(8)
            {};

        //! Suggest a new time stepsize based either on the number of newton
        //! iterations required or on the variable switch
        Scalar suggestTimeStepSize(Scalar oldTimeStep) const
            {
                /*
                  if (switched_) {
                  return somethingSmall;
                  }
                */
                // use function of the newtoncontroller
                return ParentType::suggestTimeStepSize(oldTimeStep);
            }

        //! Returns true if another iteration should be done.
        bool newtonProceed(Function &u)
            {
                return ParentType::newtonProceed(u);

                if (ParentType::newtonConverged()){
                    ParentType::method_->model().clearSwitched();
                    return false;
                }
                if (ParentType::method_->model().checkSwitched() && ParentType::numSteps_ <= minStepsSwitch && !ParentType::newtonConverged())
                    return true; // do at least some iterations after variable switch
                else if (ParentType::method_->model().checkSwitched() && ParentType::numSteps_ > minStepsSwitch) {
                    ParentType::method_->model().clearSwitched();
                    return false; // if after some iterations no convergence was reached
                }
                else
                    return ParentType::newtonProceed(u);
            }

        void newtonBeginStep()
            {
                ParentType::method_->model().setSwitchedLocalToGlobal();
            }


    protected:
        friend class NewtonControllerBase<NewtonMethod, ThisType>;
        int minStepsSwitch;

        //! called by the base class the get an indication of how physical
        //! an iterative solution is 1 means "completely physical", 0 means
        //! "completely unphysical"
        Scalar physicalness_(Function &u)
            {
                return 1.0;

                const Scalar switchVarNormFactor = 1.0; // standarization value

                // the maximum distance of a Sn value to a physically
                // meaningful value.
                Scalar maxSwitchVarDelta = 0;
                Scalar maxPwDelta = 0;
                Scalar switchVar;
                Scalar pW;

                for (int idx = 0; idx < (int) (*u).size(); idx++)
                {
                    pW = (*u)[idx][0];
                    switchVar = (*u)[idx][1];

                    if (switchVar < 0) {
                        maxSwitchVarDelta = std::max(maxSwitchVarDelta, std::abs(switchVar));
                    }
                    else if (switchVar > 1) {
                        maxSwitchVarDelta = std::max(maxSwitchVarDelta, std::abs(switchVar - 1));
                    }
                    if (pW < 0.0){
                    	maxPwDelta = std::max(maxPwDelta, std::abs(pW/1e5));
                    }
                }

                // we accept solutions up to 0.2 percent bigger than 1
                // or smaller than 0 as being physical for numerical
                // reasons...
                Scalar phys = 1.02 - maxSwitchVarDelta/switchVarNormFactor - maxPwDelta;

                // we never return exactly zero, since we want to
                // allow solutions which are "very close" to a
                // physically meaningful one
                return std::min(1.0, phys);
            }
    };
}

#endif
