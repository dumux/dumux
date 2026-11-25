// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Helper functions for normals
 */
#ifndef DUMUX_GEOMETRY_NORMAL_HH
#define DUMUX_GEOMETRY_NORMAL_HH

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Create a vector normal to the given one (v is expected to be non-zero)
 * \note This returns some orthogonal vector with arbitrary length
 * \note Components of v that are round-off compared to the largest one do not change the direction of the result
 */
template<class Vector>
inline Vector normal(const Vector& v)
{
    static_assert(Vector::size() > 1, "normal expects a coordinate dimension > 1");

    if constexpr (Vector::size() == 2)
        return Vector({-v[1], v[0]});

    using std::abs; using std::sqrt;
    using Scalar = std::decay_t<decltype(v[0])>;
    const auto maxComponent = abs(*std::max_element(v.begin(), v.end(), [](const auto& a, const auto& b){ using std::abs; return abs(a) < abs(b); }));
    const auto threshold = sqrt(std::numeric_limits<Scalar>::epsilon())*maxComponent;
    const auto it = std::find_if(v.begin(), v.end(), [&](const auto& x){ using std::abs; return abs(x) > threshold; });
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
