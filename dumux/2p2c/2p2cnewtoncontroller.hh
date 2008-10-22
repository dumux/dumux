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
        typedef TwoPTwoCNewtonController<NewtonMethod>            ThisType;
        typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

        typedef typename ParentType::Scalar            Scalar;
        typedef typename ParentType::Function          Function;
        typedef typename ParentType::JacobianAssembler JacobianAssembler;

        TwoPTwoCNewtonController(bool switched = false,
							Scalar tolerance = 1e-5,
                             int targetSteps = 8,
                             int maxSteps = 12)
            : ParentType(tolerance, targetSteps, maxSteps), switched_(switched)
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
                if (switched_ && ParentType::_numSteps < 5)
                    return true; // we always do at least two iterations

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
                const Scalar SnNormFactor = 3e-5; // reference value

                // the maximum distance of a Sn value to a physically
                // meaningful value.
                Scalar maxSnDelta = 0;
                Scalar Sn;
//                Scalar pW;

                for (int idx = 0; idx < (*u).size(); idx++)  {
//                    pW = (*u)[idx][0];
                    Sn = (*u)[idx][1];

                    if (Sn < 0) {
                        maxSnDelta = std::max(maxSnDelta, std::abs(Sn));
                    }
                    else if (Sn > 1) {
                        maxSnDelta = std::max(maxSnDelta, std::abs(Sn - 1));
                    }
//                    if (pW < 0.0){
//                    	maxSnDelta = std::max(maxSnDelta, std::abs(pW/1e5));
//                    }
                    // (so far we ignore the wetting phase pressure)
                }

                // we accept solutions up to 0.2 percent bigger than 1
                // or smaller than 0 as being physical for numerical
                // reasons...
                Scalar phys = 1.002 - maxSnDelta/SnNormFactor;

                // we never return exactly zero, since we want to
                // allow solutions which are "very close" to a
                // physically meaningful one
                return std::min(1.0, phys);
            }
        bool switched_;
    };
}

#endif
