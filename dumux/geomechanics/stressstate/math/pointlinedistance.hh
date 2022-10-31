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
 * \ingroup Geomechanics
 * \brief Stress drop law checks if the shear failure occurs and reduces
 *  the shear stress directly from the failure stress state.
 *  The failure stress state maybe rotated.
 */
#ifndef DUMUX_STRESS_STATE_MATH_POINT_LINE_DISTANCE_HH
#define DUMUX_STRESS_STATE_MATH_POINT_LINE_DISTANCE_HH

#include <cmath>

namespace Dumux{
/*!
 * \brief distance from one point to a line
 *
 * \tparam Scalar
 * \tparam Point
 * \tparam Line
 * \param p Point
 * \param l Line
 * \return Scalar
 */
template<class Scalar, class Point, class Line>
Scalar pointLineDistance(const Point& p,
                         const Line& l)
{
    using std::abs, std::sqrt;
    Scalar x = p.x();
    Scalar y = p.y();
    Scalar a = l.slope();
    Scalar b = l.intercept();
    return abs(a*x - y + b)/sqrt(1+a*a);
}
}// end namespace Dumux
#endif
