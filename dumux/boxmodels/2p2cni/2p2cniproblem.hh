// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 *
 * \brief Base class for all problems which use the non-isothermal two-phase,
 *        two-component box model
 */
#ifndef DUMUX_2P2CNI_PROBLEM_HH
#define DUMUX_2P2CNI_PROBLEM_HH

#include <dumux/boxmodels/2p2c/2p2cproblem.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCNIModel
 * \brief Base class for all problems which use the non-isothermal
 *         two-phase, two-component box model.
 */
template<class TypeTag>
class TwoPTwoCNIProblem : public TwoPTwoCProblem<TypeTag>
{
    typedef TwoPTwoCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoPTwoCNIProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
    }

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
