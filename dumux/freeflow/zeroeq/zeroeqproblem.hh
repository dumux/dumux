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
 * \brief Base class for all ZeroEq problems which use the box scheme.
 */
#ifndef DUMUX_ZEROEQ_PROBLEM_HH
#define DUMUX_ZEROEQ_PROBLEM_HH

#include <dumux/freeflow/stokes/stokesmodel.hh>
#include <dumux/freeflow/zeroeq/zeroeqproperties.hh>

namespace Dumux
{
/*!
 * \ingroup BoxZeroEqProblems
 * \brief Base class for all problems which use the ZeroEq box model.
 *
 * This implements the call of update functions used for the wall properties of
 * the ZeroEq box model.
 */
template<class TypeTag>
class ZeroEqProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    ZeroEqProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    { }

    /*!
     * \brief Called by the Dumux::TimeManager in order to initialize the problem.
     *
     * If you overload this method don't forget to call ParentType::init().<br>
     * This initializes all wall-related properties, which are necessary for the
     * ZeroEq box model.
     */
    void init()
    {
        // set the initial condition of the model
        ParentType::init();
        this->model().resetWallProperties();
        this->model().resetWallFluidProperties();
    }

    /*!
     * \brief Called by the time manager before the time integration.
     *
     * This updates all wall-related properties, which are necessary for the
     * ZeroEq box model
     */
    void preTimeStep()
    {
        this->model().updateWallProperties();
    }
};

}

#endif // DUMUX_ZEROEQ_PROBLEM_HH
