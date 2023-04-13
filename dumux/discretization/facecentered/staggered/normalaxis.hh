// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_NORMAL_AXIS_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_NORMAL_AXIS_HH

#include <cstddef>
#include <algorithm>
#include <numeric>

namespace Dumux {

/*!
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
