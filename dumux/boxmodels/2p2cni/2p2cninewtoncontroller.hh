// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
 * \brief A 2p2cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P2CNI_NEWTON_CONTROLLER_HH
#define DUMUX_2P2CNI_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup TwoPTwoCNIModel
 * \brief A 2p2cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPTwoCNINewtonController : public NewtonController<TypeTag >
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;


    typedef typename SolutionTypes::SolutionFunction SolutionFunction;

public:
    TwoPTwoCNINewtonController()
    {
        this->setRelTolerance(1e-5);
        this->setTargetSteps(9);
        this->setMaxSteps(18);
    };

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    void newtonEndStep(SolutionFunction &u, SolutionFunction &uOld)
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

    //! Returns true if another iteration should be done.
    bool newtonProceed(SolutionFunction &u)
    {
        return ParentType::newtonProceed(u);

        bool baseProceed = ParentType::newtonProceed(u);

        // if we just switched some primary variables we
        // proceed for at least for newton iterations if
        // before we give up.
        if (!baseProceed &&
            ParentType::model().switched() &&
            ParentType::numSteps_ < 4 &&
            !ParentType::newtonConverged())
        {
            return true;
        }

        return baseProceed;
    }

protected:
    friend class NewtonController<TypeTag>;
    //! called by the base class to get an indication of how physical
    //! an iterative solution is 1 means "completely physical", 0 means
    //! "completely unphysical"
    Scalar physicalness_(SolutionFunction &u)
    {
        return 1.0;

        const Scalar switchVarNormFactor = 1e-1; // standardisation value

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

        // we accept solutions up to 1 percent bigger than 1
        // or smaller than 0 as being physical for numerical
        // reasons...
        Scalar phys = 1.01 - maxSwitchVarDelta/switchVarNormFactor - maxPwDelta;

        return std::min(1.0, phys);
    }
};
}

#endif
