// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 *
 * \brief Base class for all problems which use the non-isothermal
 *        two-phase box model
 */
#ifndef DUMUX_2PNI_PROBLEM_HH
#define DUMUX_2PNI_PROBLEM_HH

#include <dumux/boxmodels/2p/2pproblem.hh>


namespace Dumux
{
/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxBaseProblems
 * \brief Base class for all problems which use the non-isothermal
 *         two-phase box model.
 *
 */
template<class TypeTag>
class TwoPNIProblem : public TwoPProblem<TypeTag>
{
    typedef TwoPProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
        */
    TwoPNIProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This is a non-isothermal model, so this method just throws an
     * exception. This method MUST NOT be overwritten by the actual
     * problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::Exception, "temperature() method called for a 2p2cni problem"); };

    // \}
};

}

#endif
