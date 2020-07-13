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
 * \brief This file contains functions related to calculate pore body properties.
 *
 */
#ifndef DUMUX_PNM_BASE_PORE_PROPERTIES_HH
#define DUMUX_PNM_BASE_PORE_PROPERTIES_HH

#include <string>

#include <dune/common/exceptions.hh>

namespace Dumux {

namespace Pore {

enum class Shape
{ circle, square, cube, sphere, cylinder };


//! Get the shape from a string description of the shape
inline std::string shapeToString(Shape s)
{
    switch (s)
    {
        case Shape::square: return "Square";
        case Shape::circle: return "Circle";
        case Shape::cube: return "Cube";
        case Shape::sphere: return "Sphere";
        case Shape::cylinder: return "Cylinder";
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape!");
    }
}

//! Get the shape from a string description of the shape
inline Shape shapeFromString(const std::string& s)
{
    if (s == shapeToString(Shape::square)) return Shape::square;
    else if (s == shapeToString(Shape::circle)) return Shape::circle;
    else if (s == shapeToString(Shape::cube)) return Shape::cube;
    else if (s == shapeToString(Shape::sphere)) return Shape::sphere;
    else if (s == shapeToString(Shape::cylinder)) return Shape::cylinder;
    else DUNE_THROW(Dune::InvalidStateException, s << " is not a valid shape");
}

/*!
* \brief Returns the volume of a given geometry based on the inscribed radius
*/
template<class Scalar>
inline Scalar volume(Shape shape, Scalar inscribedRadius)
{
    switch(shape)
    {
        case Shape::cube: return 8*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::sphere: return 4.0/3.0*M_PI*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::circle: return M_PI*inscribedRadius*inscribedRadius; break;
        case Shape::square: return 4.0*inscribedRadius*inscribedRadius; break;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry");
    }
}

/*!
* \brief Returns the volume of a given geometry based on the inscribed radius and the height
*/
template<class Scalar>
inline Scalar volume(Shape shape, Scalar inscribedRadius, Scalar height)
{
    switch(shape)
    {
        case Shape::cylinder: return M_PI*inscribedRadius*inscribedRadius*height; break;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry");
    }
}

} // end namespace
}

#endif
