// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Compute integrals over serveral geometry types
 */
#ifndef DUMUX_GEOMETRY_INTEGRATE_HH
#define DUMUX_GEOMETRY_INTEGRATE_HH

#include <cmath>
#include <limits>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Integrate a function over some geometry
 */
template<class Function, class Geometry>
auto integrate(const Function& f, const Geometry& geometry, unsigned int integrationOrder = 4)
{
    using Scalar = typename Geometry::ctype;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using FunctionType = std::decay_t<decltype(f(GlobalPosition(0.0)))>;
    FunctionType integral(0.0);
    const auto &quad = Dune::QuadratureRules<Scalar, Geometry::mydimension>::rule(geometry.type(), integrationOrder);
    for (auto &&qp : quad)
    {
        const auto ipGlobal = geometry.global(qp.position());
        integral += qp.weight() * geometry.integrationElement(qp.position()) * f(ipGlobal);
    }

    return integral;
}

/*!
 * \ingroup Geometry
 * \brief Average integral of a function over some geometry
 */
template<class Function, class Geometry>
auto integralAverage(const Function& f, const Geometry& geometry, unsigned int integrationOrder = 4)
{
    using Scalar = typename Geometry::ctype;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using FunctionType = std::decay_t<decltype(f(GlobalPosition(0.0)))>;
    FunctionType integral(0.0);
    Scalar volume(0.0);
    const auto &quad = Dune::QuadratureRules<Scalar, Geometry::mydimension>::rule(geometry.type(), integrationOrder);
    for (auto &&qp : quad)
    {
        const auto ipGlobal = geometry.global(qp.position());
        Scalar vol = qp.weight() * geometry.integrationElement(qp.position());
        integral += vol * f(ipGlobal);
        volume += vol;
    }
    integral /= volume;

    return integral;
}

} // end namespace Dumux

#endif
