// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief This file contains functions related to calculate pore-throat properties.
 */
#ifndef DUMUX_PNM_BASE_THROAT_PROPERTIES_HH
#define DUMUX_PNM_BASE_THROAT_PROPERTIES_HH

#include <string>
#include <cmath>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

namespace Dumux::PoreNetwork::Throat {

//! Collection of different pore-throat shapes
enum class Shape
{ scaleneTriangle, equilateralTriangle, square, rectangle, circle, twoPlates, polygon };

//! Get the shape from a string description of the shape
inline std::string shapeToString(Shape s)
{
    switch (s)
    {
        case Shape::scaleneTriangle: return "ScaleneTriangle";
        case Shape::equilateralTriangle: return "EquilateralTriangle";
        case Shape::square: return "Square";
        case Shape::rectangle: return "Rectangle";
        case Shape::circle: return "Circle";
        case Shape::twoPlates: return "TwoPlates";
        case Shape::polygon: return "Polygon";
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape!");
    }
}

//! Get the shape from a string description of the shape
inline Shape shapeFromString(const std::string& s)
{
    if (s == shapeToString(Shape::equilateralTriangle)) return Shape::equilateralTriangle;
    else if (s == shapeToString(Shape::scaleneTriangle)) return Shape::scaleneTriangle;
    else if (s == shapeToString(Shape::square)) return Shape::square;
    else if (s == shapeToString(Shape::rectangle)) return Shape::rectangle;
    else if (s == shapeToString(Shape::circle)) return Shape::circle;
    else if (s == shapeToString(Shape::twoPlates)) return Shape::twoPlates;
    else if (s == shapeToString(Shape::polygon)) return Shape::polygon;
    else DUNE_THROW(Dune::InvalidStateException, s << " is not a valid shape");
}

/*!
 * \brief Returns the radius of a pore throat
 *
 * \param poreRadiusOne The radius of the first pore
 * \param poreRadiusTwo The radius of the second pore
 * \param centerTocenterDist The center-to-center distance between the pores
 * \param n Fitting parameter
 *
 * Joekar-Niasar et al. (2008) https://doi.org/10.1007/s11242-007-9191-7
 */
template<class Scalar>
inline Scalar averagedRadius(const Scalar poreRadiusOne, const Scalar poreRadiusTwo, const Scalar centerTocenterDist, const Scalar n = 0.1)
{
   assert(n > 0.0);
   using std::sin; using std::cos; using std::pow;
   const Scalar rOneTilde = poreRadiusOne/centerTocenterDist;
   const Scalar rTwoTilde = poreRadiusTwo/centerTocenterDist;
   const Scalar a = sin(M_PI/4.0);
   const Scalar b = cos(M_PI/4.0);
   const Scalar rhoOne = rOneTilde*a / pow((1.0 - rOneTilde*b), n);
   const Scalar rhoTwo = rTwoTilde*a / pow((1.0 - rTwoTilde*b), n);
   const Scalar rTilde = rhoOne*rhoTwo * pow((pow(rhoOne, 1.0/n) + pow(rhoTwo, 1.0/n)), -n);
   return rTilde * centerTocenterDist;
}


//! Returns the corner half angle
template<class Scalar>
inline Dune::ReservedVector<Scalar, 4> cornerHalfAngles(Shape shape)
{
    switch(shape)
    {
        case Shape::equilateralTriangle:
        {
            const Scalar value = M_PI/6.0;
            return Dune::ReservedVector<Scalar, 4>{value, value, value};
        }
        case Shape::square:
        {
            const Scalar value = M_PI/4.0;
            return Dune::ReservedVector<Scalar, 4>{value, value, value, value};
        }
        case Shape::rectangle:
        {
            const Scalar value = M_PI/4.0;
            return Dune::ReservedVector<Scalar, 4>{value, value, value, value};
        }
        case Shape::circle:
        {
            const Scalar value = 0.5*M_PI; // we define the (single) corner angle of a circle as 180°
            return { value };
        }
        case Shape::polygon: DUNE_THROW(Dune::InvalidStateException, "Corner half angles for polygons must be calculated explicitly");
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape");
        // TODO implement angles for scaleneTriangle as given by Valvatne & Blunt (2004)
    }
}

//! Returns the value of the shape factor for an equilateral triangle
template<class Scalar>
inline constexpr Scalar shapeFactorEquilateralTriangle() noexcept
{
    using std::sqrt;
    return sqrt(3.0)/36.0;
}

//! Returns the value of the shape factor for a square
template<class Scalar>
inline constexpr Scalar shapeFactorSquare() noexcept
{
    return 1.0/16.0;
}

//! Returns the value of the shape factor for a rectangle
template<class Scalar>
inline constexpr Scalar shapeFactorRectangle(const Scalar inscribedRadius, const Scalar height) noexcept
{
    const Scalar a = 2.0*inscribedRadius; // shorter side length
    const Scalar b = height; // longer side length
    const Scalar area = a*b;
    const Scalar perimeter = 2.0*a + 2.0*b;
    return area / (perimeter*perimeter);
}

//! Returns the value of the shape factor for a circle
template<class Scalar>
inline constexpr Scalar shapeFactorCircle() noexcept
{
    return 1.0/(4.0*M_PI);
}

//! Returns the value of the shape factor for a given shape
template<class Scalar>
inline Scalar shapeFactor(Shape shape, const Scalar inscribedRadius)
{
    switch(shape)
    {
        case Shape::equilateralTriangle: return shapeFactorEquilateralTriangle<Scalar>();
        case Shape::square: return shapeFactorSquare<Scalar>();
        case Shape::circle: return shapeFactorCircle<Scalar>();
        case Shape::twoPlates: return 0.0; // TODO is this a good idea?
        case Shape::polygon: DUNE_THROW(Dune::InvalidStateException, "The shape factor for a polygon has to be defined by the input data");
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown shape");
    }
}

//! Returns the shape for a given shape factor
template<class Scalar>
inline constexpr Shape shape(const Scalar shapeFactor) noexcept
{
    if (shapeFactor < shapeFactorEquilateralTriangle<Scalar>())
        return Shape::scaleneTriangle;
    else if (shapeFactor <= shapeFactorEquilateralTriangle<Scalar>())
        return Shape::equilateralTriangle;
    else if (shapeFactor < shapeFactorSquare<Scalar>())
        return Shape::rectangle;
    else if (shapeFactor == shapeFactorSquare<Scalar>())
        return Shape::square;
    else if (shapeFactor == shapeFactorCircle<Scalar>())
        return Shape::circle;
    else
        return Shape::polygon;
}

//! Returns if a shape is regular
inline bool isRegularShape(Shape shape)
{
    switch (shape)
    {
        case Shape::equilateralTriangle: return true;
        case Shape::square: return true;
        case Shape::circle: return true;
        case Shape::rectangle: return false;
        case Shape::scaleneTriangle: return false;
        case Shape::polygon: DUNE_THROW(Dune::InvalidStateException, "Equality of Corner half angles for polygons must be determined explicitly");
        default:
            DUNE_THROW(Dune::InvalidStateException, "Unknown shape");
    }
}

//! Returns the cross-sectional area of a given geometry
template<class Scalar>
inline Scalar totalCrossSectionalArea(const Shape shape, const Scalar inscribedRadius)
{
    using std::sqrt;
    switch(shape)
    {
        case Shape::equilateralTriangle: return 3.0*sqrt(3.0)*inscribedRadius*inscribedRadius;
        case Shape::square: return 4.0*inscribedRadius*inscribedRadius;
        case Shape::circle: return M_PI*inscribedRadius*inscribedRadius;
        case Shape::twoPlates: return 2.0*inscribedRadius;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry: " << shapeToString(shape));
    }
}

//! Returns the cross-sectional area of a rectangle
template<class Scalar>
inline constexpr Scalar totalCrossSectionalAreaForRectangle(const Scalar inscribedRadius, const Scalar height) noexcept
{
    return 2.0*inscribedRadius * height;
}

//! Returns the number of corners of a given geometry
inline std::size_t numCorners(Shape shape)
{
    switch(shape)
    {
        case Shape::scaleneTriangle: return 3;
        case Shape::equilateralTriangle: return 3;
        case Shape::square: return 4;
        case Shape::rectangle: return 4;
        case Shape::circle: return 0;
        default : DUNE_THROW(Dune::InvalidStateException, "Unsupported geometry: " << shapeToString(shape));
    }
}

/*!
 * \brief Return the cross-sectional area of a wetting layer residing in a corner of a throat
 *
 * \param curvatureRadius The radius of curvature within the throat (surfaceTension/pc)
 * \param contactAngle The contact angle within the throat
 * \param cornerHalfAngle The corner half-angle of the the throat's corner of interest
 *
 * See, e.g, Valvatne & Blunt (2004), eq. A7 https://doi.org/10.1029/2003WR002627
 */
template<class Scalar>
inline constexpr Scalar wettingLayerCrossSectionalArea(const Scalar curvatureRadius,
                                                       const Scalar contactAngle,
                                                       const Scalar cornerHalfAngle) noexcept
{
    using std::sin;
    using std::cos;
    return curvatureRadius*curvatureRadius *(cos(contactAngle) * cos(contactAngle + cornerHalfAngle) / sin(cornerHalfAngle)
           + cornerHalfAngle + contactAngle - M_PI/2.0);
}

} // end Dumux::PoreNetwork::Throat

#endif
