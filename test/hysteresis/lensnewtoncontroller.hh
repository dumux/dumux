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
 * \brief A Pw-Sn specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_LENS_NEWTON_CONTROLLER_HH
#define DUMUX_LENS_NEWTON_CONTROLLER_HH

#include <dumux/new_models/box/pwsn/pwsnnewtoncontroller.hh>

namespace Dune {
namespace Lens {
    /*!
     * \brief A Pw-Sn specific controller for the newton solver for
     * the lens problem.
     *
     * This  controller  is   basically  a  PwSnNewtonController,  but
     * notifies the LensSimulation class  when a newton step has ended
     * so that  the simulation has the  chance to write  some stuff to
     * disc which aids to analyze the convergence behaviour.
     */
    template <class NewtonMethod, class Simulation>
    class LensNewtonController
        : public PwSnNewtonController<NewtonMethod>
    {
    public:
        typedef PwSnNewtonController<NewtonMethod>     ParentType;

        typedef typename ParentType::Scalar            Scalar;
        typedef typename ParentType::Function          Function;
        typedef typename ParentType::JacobianAssembler JacobianAssembler;

        LensNewtonController(Simulation &sim,
                             Scalar tolerance = 1e-5,
                             int targetSteps = 6,
                             int maxSteps = 12)
            : ParentType(tolerance, targetSteps, maxSteps), _sim(sim)
            {};

        //! Indicates that the newton method is started.
        void newtonBegin(NewtonMethod *method, Function &u)
            {
                // notify the PwSnNewtonController
                ParentType::newtonBegin(method, u);

                // notify the simulation controller
                _sim.newtonBegin();
            }

        //! Indicates that we're done solving one newton step.
        void newtonEndStep(Function &u, Function &uOld)
            {
                // notify the PwSnNewtonController
                ParentType::newtonEndStep(u, uOld);

                // notify the simulation controller
                _sim.newtonEndStep(u, uOld);
            };

        //! Indicates that the newton method has finished.
        void newtonEnd()
            {
                // notify the PwSnNewtonController
                ParentType::newtonEnd();

                // notify the simulation controller
                _sim.newtonEnd();
            }

    private:
        Simulation &_sim;
    };
}
}

#endif
