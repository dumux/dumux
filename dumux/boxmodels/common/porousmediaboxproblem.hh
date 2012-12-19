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
 * \brief Base class for all problems which use the two-phase box model
 */
#ifndef DUMUX_POROUS_MEDIA_BOX_PROBLEM_HH
#define DUMUX_POROUS_MEDIA_BOX_PROBLEM_HH

#include <dumux/implicit/common/implicitporousmediaproblem.hh>

namespace Dumux
{
/*!
 * \ingroup BoxBaseProblems
 * \brief Base class for all porous media box problems
 */
template<class TypeTag>
class PorousMediaBoxProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    DUNE_DEPRECATED_MSG("Use ImplicitPorousMediaProblem from dumux/implicit/common/implicitporousmediaproblem.hh.")
    PorousMediaBoxProblem(TimeManager &timeManager,
                const GridView &gridView,
                const bool verbose = true)
        : ParentType(timeManager, gridView)
    {}
};

}

#endif
