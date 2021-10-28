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
 * \brief Define geometric tolerances.
 */
#ifndef DUMUX_GEOMETRY_PRECISION_HH
#define DUMUX_GEOMETRY_PRECISION_HH

#include <string>
#include <dumux/common/parameters.hh>

namespace Dumux::Geometry {

/*!
 * \ingroup Geometry
 * \brief Define geometric tolerances.
 */
template<class ctype = double>
class Precision
{
public:
    /*!
     * \brief Return the relative tolerance.
     */
    static ctype relativeTolerance()
    {
        static const ctype relTolerance = getParam<ctype>("Geometry.RelativeTolerance");
        return relTolerance;
    }
};

} // end namespace Dumux::Geometry

#endif
