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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief A function to compute bounding spheres of points clouds or convex polytopes
 */
#ifndef DUMUX_GEOMETRY_SPHERE_HH
#define DUMUX_GEOMETRY_SPHERE_HH

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A simple sphere type
 */
template<class Scalar, int dim>
class Sphere
{
public:
    using Point = Dune::FieldVector<Scalar, dim>;

    Sphere()
    : center_(0.0)
    , radius_(0.0)
    {}

    Sphere(const Point& center, Scalar radius)
    : center_(center)
    , radius_(radius)
    {}

    Scalar radius() const
    { return radius_; }

    const Point& center() const
    { return center_; }

private:
    Point center_;
    Scalar radius_;
};

} // end namespace Dumux

#endif
