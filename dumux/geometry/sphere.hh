// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief A function to compute bounding spheres of points clouds or convex polytopes
 */
#ifndef DUMUX_GEOMETRY_SPHERE_HH
#define DUMUX_GEOMETRY_SPHERE_HH

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A simple sphere type
 */
template<class Scalar, int dim>
class Sphere
{
public:
    using Point = Dune::FieldVector<Scalar, dim>;

    Sphere()
    : center_(0.0)
    , radius_(0.0)
    {}

    Sphere(const Point& center, Scalar radius)
    : center_(center)
    , radius_(radius)
    {}

    Scalar radius() const
    { return radius_; }

    const Point& center() const
    { return center_; }

private:
    Point center_;
    Scalar radius_;
};

} // end namespace Dumux

#endif
