//$Id: 2p2cninewtoncontroller.hh 768 2008-11-05 13:14:44Z melanie $
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
#ifndef DUMUX_2P2CNI_NEWTON_CONTROLLER_HH
#define DUMUX_2P2CNI_NEWTON_CONTROLLER_HH

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
        typedef TwoPTwoCNINewtonController<NewtonMethod>		ThisType;
        typedef NewtonControllerBase<NewtonMethod, ThisType>  	ParentType;
        typedef typename NewtonMethod::Model  					Model_TEst;
        typedef typename ParentType::Scalar            			Scalar;
        typedef typename ParentType::Function          			Function;
        typedef typename ParentType::JacobianAssembler 			JacobianAssembler;

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
			/*
				if (variableSwitch) {
					return somethingSmall;
				}
				*/
				// use function of the newtoncontroller
				return ParentType::suggestTimeStepSize(oldTimeStep);
			}

        //! Returns true iff another iteration should be done.
        bool newtonProceed(Function &u)
            {
//				Model_TEst *test;
//				switched_ = test->localJacobian.switched;
//				if (switched_)
//				{
//					test->localJacobian.switchBreak = true;
//
//					if(switchCount_ == 5)
//					{
//					 test->localJacobian.switchBreak = false;
//					 switchCount_ = 0;
//					 return ParentType::newtonProceed(u);
//					}
//					else
//					{
//						switchCount_++;
//						return true; // we always do at least five iterations
//					}
//				}
                return ParentType::newtonProceed(u);
            }


    protected:
        friend class NewtonControllerBase<NewtonMethod, ThisType>;
        //! called by the base class the get an indication of how physical
        //! an iterative solution is 1 means "completely physical", 0 means
        //! "completely unphysical"
        Scalar _physicalness(Function &u)
            {
//                return 1.0;
                const Scalar SnNormFactor = 1e-5; // reference value
                const Scalar pWNormFactor = 1.;
                const Scalar tempNormFactor = 1e-3;
                // the maximum distance of a Sn value to a physically
                // meaningful value.
                Scalar maxSwitchVarDelta = 0;
                Scalar maxPWDelta = 0;
                Scalar maxTempDelta = 0;
                Scalar switchVar;
                Scalar pW;
                Scalar temp;
//                Scalar pW;

                for (int idx = 0; idx < (*u).size(); idx++)  {

                	pW = (*u)[idx][0];
                    switchVar = (*u)[idx][1];
                    temp = (*u)[idx][2];

                    if (pW < 1.e5) {
                    	maxPWDelta = std::max(maxPWDelta, std::abs(pW-1.e5));
                    }
                    if (pW > 2.e8) {
                    	maxPWDelta = std::max(maxPWDelta, std::abs(pW-2.e8));
                    }
                    if (switchVar < 0) {
                        maxSwitchVarDelta = std::max(maxSwitchVarDelta, std::abs(switchVar));
                    }
                    if (switchVar > 1) {
                        maxSwitchVarDelta = std::max(maxSwitchVarDelta, std::abs(switchVar - 1));
                    }
                    if (temp < 273.) {
                    	maxTempDelta = std::max(maxTempDelta, std::abs(temp-273.));
                    }
                    if (temp > 600.) {
                    	maxTempDelta = std::max(maxTempDelta, std::abs(temp-600.));
                    }


                }

                // we accept solutions up to 0.2 percent bigger than 1
                // or smaller than 0 as being physical for numerical
                // reasons...

                Scalar phys = 1.002 - maxSwitchVarDelta/SnNormFactor - maxPWDelta/pWNormFactor - maxTempDelta/tempNormFactor;

                // we never return exactly zero, since we want to
                // allow solutions which are "very close" to a
                // physically meaningful one

                return std::min(1.0, phys);
            }
        bool switched_;
        int switchCount_;
    };
}

#endif
