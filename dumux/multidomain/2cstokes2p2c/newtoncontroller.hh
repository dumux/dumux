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
 * \brief Reference implementation of a Newton controller for the coupling of a two-component Stokes model
 *        and a two-phase two-component porous-medium model under isothermal conditions.
 */
#ifndef DUMUX_2CSTOKES_2P2C_NEWTON_CONTROLLER_HH
#define DUMUX_2CSTOKES_2P2C_NEWTON_CONTROLLER_HH

#include <dumux/multidomain/newtoncontroller.hh>

namespace Dumux
{

/*!
 * \ingroup Newton
 * \ingroup TwoPTwoCStokesTwoCModel
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief Implementation of a Newton controller for the coupling of a two-component Stokes model
 *        and a two-phase two-component porous-medium model under isothermal conditions.
 *
 * The Newton controller ensures that the updateStaticData routine is called
 * in the porous-medium sub-problem
 */
template <class TypeTag>
class TwoCStokesTwoPTwoCNewtonController : public MultiDomainNewtonController<TypeTag>
{
    typedef MultiDomainNewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

public:
    //! \brief The constructor
    TwoCStokesTwoPTwoCNewtonController(const Problem &problem)
        : ParentType(problem)
    {  }

    //! \copydoc Dumux::NewtonController::newtonEndStep()
    void newtonEndStep(SolutionVector &uCurrentIter, SolutionVector &uLastIter)
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        this->model_().sdModel2().updateStaticData(this->model_().sdModel2().curSol(),
                                                   this->model_().sdModel2().prevSol());
    }
};

} // namespace Dumux

#endif // DUMUX_2CSTOKES_2P2C_NEWTON_CONTROLLER_HH
