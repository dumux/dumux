// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief A 2p1c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P1CNI_NEWTON_CONTROLLER_HH
#define DUMUX_2P1CNI_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/porousmediumflow/2p2c/implicit/newtoncontroller.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup TwoPOneCModel
 * \brief A 2p1cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPOneCNINewtonController : public TwoPTwoCNewtonController<TypeTag>
{
    typedef TwoPTwoCNewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

public:
    TwoPOneCNINewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    void newtonBeginStep()
    {
        // call the method of the base class
        ParentType::newtonBeginStep();

        //reset the indicator of spurious cold water fluxes into the steam zone
        this->model_().localJacobian().localResidual().resetSpuriousFlowDetected();
    }


    /*!
     * \brief Called after each Newton update
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     */
    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        // call the method of the base class
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        //Output message after each Newton Step if spurious fluxes have been blocked
        if (this->model_().localJacobian().localResidual().spuriousFlowDetected())
            std::cout <<"Spurious fluxes blocked!" <<std::endl;
    }

};
}

#endif
