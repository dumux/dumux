// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief This file contains functions related to calculate pore-body properties.
 */
#ifndef DUMUX_PNM_BASE_PORE_PROPERTIES_HH
#define DUMUX_PNM_BASE_PORE_PROPERTIES_HH

#include <cmath>
#include <string>
#include <dune/common/exceptions.hh>

namespace Dumux::PoreNetwork::Pore {

//! Collection of different pore-body shapes
enum class Shape
{ circle, square, cube, sphere, cylinder, tetrahedron, octahedron, icosahedron, dodecahedron };

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
        case Shape::tetrahedron: return "Tetrahedron";
        case Shape::octahedron: return "Octahedron";
        case Shape::icosahedron: return "Icosahedron";
        case Shape::dodecahedron: return "Dodecahedron";
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
    else if (s == shapeToString(Shape::tetrahedron)) return Shape::tetrahedron;
    else if (s == shapeToString(Shape::octahedron)) return Shape::octahedron;
    else if (s == shapeToString(Shape::icosahedron)) return Shape::icosahedron;
    else if (s == shapeToString(Shape::dodecahedron)) return Shape::dodecahedron;
    else DUNE_THROW(Dune::InvalidStateException, s << " is not a valid shape");
}


//! Returns the volume of a given geometry based on the inscribed radius
template<class Scalar>
inline Scalar volume(Shape shape, Scalar inscribedRadius)
{
    switch(shape)
    {
        case Shape::cube: return 8*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::sphere: return 4.0/3.0*M_PI*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::circle: return M_PI*inscribedRadius*inscribedRadius; break;
        case Shape::square: return 4.0*inscribedRadius*inscribedRadius; break;
        case Shape::tetrahedron: return 13.85*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::octahedron: return 6.93*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::icosahedron: return 5.05*inscribedRadius*inscribedRadius*inscribedRadius; break;
        case Shape::dodecahedron: return 5.55*inscribedRadius*inscribedRadius*inscribedRadius; break;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry");
    }
}


//! Returns the volume of a given geometry based on the inscribed radius and the height
template<class Scalar>
inline Scalar volume(Shape shape, Scalar inscribedRadius, Scalar height)
{
    switch(shape)
    {
        case Shape::cylinder: return M_PI*inscribedRadius*inscribedRadius*height; break;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry");
    }
}

} // end Dumux::PoreNetwork::Pore

#endif
