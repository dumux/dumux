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
 * \brief Helper functions for normals
 */
#ifndef DUMUX_GEOMETRY_NORMAL_HH
#define DUMUX_GEOMETRY_NORMAL_HH

#include <algorithm>
#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Create a vector normal to the given one (v is expected to be non-zero)
 * \note This returns some orthogonal vector with arbitrary length
 */
template<class Vector>
inline Vector normal(const Vector& v)
{
    static_assert(Vector::size() > 1, "normal expects a coordinate dimension > 1");

    if constexpr (Vector::size() == 2)
        return Vector({-v[1], v[0]});

    const auto it = std::find_if(v.begin(), v.end(), [](const auto& x) { return Dune::FloatCmp::ne(x, 0.0); });
    const auto index = std::distance(v.begin(), it);
    if (index != Vector::size()-1)
    {
        Vector normal(0.0);
        normal[index] = -v[index+1];
        normal[index+1] = v[index];
        return normal;
    }
    else
    {
        Vector normal(0.0);
        normal[index-1] = -v[index];
        normal[index] = v[index-1];
        return normal;

    }
}

/*!
 * \ingroup Geometry
 * \brief Create a vector normal to the given one (v is expected to be non-zero)
 * \note This returns some orthogonal vector with unit length
 */
template<class Vector>
inline Vector unitNormal(const Vector& v)
{
    auto normal = Dumux::normal(v);
    normal /= normal.two_norm();
    return normal;
}

} // end namespace Dumux

#endif
