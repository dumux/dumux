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
 * \brief Line class in the Mohr Space
 */
#ifndef DUMUX_STRESS_STATE_MATH_LINE_HH
#define DUMUX_STRESS_STATE_MATH_LINE_HH

namespace Dumux{
/*!
 * \brief Line class in Mohr Space
 *        in form y = ax + b,
 *        where a is slope and b is the intercept.
 *
 * \tparam Scalar
 */
template<class Scalar>
class Line{
public:
    Line(const Scalar& slope, const Scalar& intercept)
    :slope_(slope),intercept_(intercept)
    {}

    Line():Line(0.0,0.0){}

    Scalar slope() const { return slope_; }
    Scalar intercept() const { return intercept_; }
    Scalar y(const Scalar& x) const { return slope_ * x + intercept_; }
private:
    Scalar slope_;
    Scalar intercept_;
};
}// end namespace Dumux
#endif
