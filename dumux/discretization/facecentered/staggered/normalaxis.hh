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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_NORMAL_AXIS_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_NORMAL_AXIS_HH

#include <cstddef>
#include <algorithm>

namespace Dumux {

/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Returns the normal axis index of a unit vector (0 = x, 1 = y, 2 = z)
 */
template<class Vector>
inline static std::size_t normalAxis(const Vector& v)
{
    using std::abs;

    constexpr auto eps = 1e-8;
    assert(std::any_of(v.begin(), v.end(), [=](auto x){ return abs(x) > eps; }));

    const auto result = std::distance(
        std::begin(v), std::find_if(v.begin(), v.end(), [eps=eps](auto x){ return abs(x) > eps; })
    );

    // make sure there is only one non-zero entry
    assert(v[result] == std::accumulate(v.begin(), v.end(), 0.0));

    return result;
}

} // end namespace Dumux

#endif
