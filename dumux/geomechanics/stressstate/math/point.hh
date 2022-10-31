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
 * \brief Point class in Mohr Space
 */
#ifndef DUMUX_STRESS_STATE_MATH_POINT_HH
#define DUMUX_STRESS_STATE_MATH_POINT_HH

#include <dune/common/fvector.hh>

namespace Dumux{
template<class Scalar>
class Point{
    using Position = Dune::FieldVector<Scalar,2>;
public:
    Point(const Scalar& x, const Scalar& y)
    {
        pos_[0] = x;
        pos_[1] = y;
    }

    Point():Point(0.0,0.0){};
    Scalar x() const { return pos_[0]; }
    Scalar y() const { return pos_[1]; }

    /*!
     * \brief return the position as vector
     *
     * \return Position
     */
    Position pos() const {return pos_;}
private:
    Position pos_;
};
}// end namespace Dumux
#endif
