// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
