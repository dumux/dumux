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
 * \brief This file contains functions useful for all types of pore-network models.
 *
 */
#ifndef DUMUX_PNM_GEOMETRY_HH
#define DUMUX_PNM_GEOMETRY_HH

#include <cmath>
#include <dune/common/exceptions.hh>


namespace Dumux {

/*!
* \brief Returns the projected radius of a throat cutting a plane
*
* \param throatRadius The throat's actual radius
* \param element The element
* \param cutPlaneNormal The plane's normal vector
*
*/
template<class Scalar, class Element>
inline Scalar projectedThroatRadius(const Scalar& throatRadius,
                                    const Element& element,
                                    const typename Element::Geometry::GlobalCoordinate& cutPlaneNormal)
{
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    if (GlobalPosition::dimension != 2)
        DUNE_THROW(Dune::InvalidStateException, "Projection of throat radius only works in 2D");

    // determine the throat radius, which might depend on the throat's angle of orientation
    const GlobalPosition throatOrientationVector = element.geometry().corner(1) - element.geometry().corner(0);

    // get the angle between the throat and the plane it cuts
    using std::acos;
    using std::sin;
    const Scalar cosPhi = cutPlaneNormal * throatOrientationVector / (cutPlaneNormal.two_norm() * throatOrientationVector.two_norm());
    const Scalar phi = acos(cosPhi); // in radians
    const Scalar alpha = phi > 0.5*M_PI ? phi - 0.5*M_PI : 0.5*M_PI - phi; // 0.5*pi rad == 90 deg; we define alpha always smaller than 90 deg

    // get the radius projected onto the plane the throat cuts
    return throatRadius / sin(alpha);
}

} // end namespace

#endif
