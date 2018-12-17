// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup CCDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_CC_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_CC_ELEMENT_BOUNDARY_TYPES_HH

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Boundary types gathered on an element
 * \note This class exists only for compatibility purposes with the
 *        box scheme. The cell-centered schemes and the box scheme use
 *        a common base local residual, which passes an ElementBoundaryTypes
 *        object to the implemented interfaces.
 */
class CCElementBoundaryTypes
{
public:
    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element finite volume geometry
     */
    template<class Problem, class Element, class FVElementGeometry>
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {}
};

} // namespace Dumux

#endif
