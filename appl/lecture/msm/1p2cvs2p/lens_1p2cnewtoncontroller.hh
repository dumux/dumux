// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief A 1p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_LENS_1P2C_NEWTON_CONTROLLER_HH
#define DUMUX_LENS_1P2C_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup TwoPTwoCBoxModel
 * \brief A 2p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class LensOnePTwoCNewtonController
    : public NewtonController<TypeTag >
{
    typedef LensOnePTwoCNewtonController<TypeTag>  ThisType;
    typedef NewtonController<TypeTag>      ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;

public:
    LensOnePTwoCNewtonController()
    {
        this->setRelTolerance(1e-7);
        this->setTargetSteps(9);
        this->setMaxSteps(18);

        //load interface-file

        Dumux::InterfaceProblemProperties interfaceProbProps("interface1p2c.xml");

        maxTimeStepSize_ = interfaceProbProps.IPP_MaxTimeStepSize;
    };

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // use function of the newtoncontroller
        return std::min(maxTimeStepSize_, ParentType::suggestTimeStepSize(oldTimeStep));
    }

private:
    Scalar maxTimeStepSize_;
};
}

#endif
