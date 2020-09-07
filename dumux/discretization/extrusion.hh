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
 * \ingroup Discretization
 * \brief Helper classes to compute the integration elements
 *
 * Provides area, volume, and integration elements for integration formulas.
 * Modifying these quantities is useful to realize extrusions of the computational domain.
 */
#ifndef DUMUX_DISCRETIZATION_EXTRUSION_HH
#define DUMUX_DISCRETIZATION_EXTRUSION_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Default implementation that performs no extrusion (extrusion with identity)
 */
struct NoExtrusion
{
    template<class SCVF>
    static constexpr auto area(const SCVF& scvf)
    { return scvf.area(); }

    template<class SCV>
    static constexpr auto volume(const SCV& scv)
    { return scv.volume(); }

    template<class Geometry>
    static constexpr auto integrationElement(const Geometry& geo, const typename Geometry::LocalCoordinate& x)
    { return geo.integrationElement(x); }
};

/*!
 * \ingroup Discretization
 * \brief Rotation symmetric extrusion policy for rotating about an external axis
 * \tparam radAx The radial axis perpendicular to the symmetry axis (0 = x, 1 = y)
 */
template<int radAx = 0>
struct RotationalExtrusion
{
    static constexpr int radialAxis = radAx;

    /*!
     * \brief Transformed sub-control-volume face area
     * \note Mid-point rule integrals are only exact for constants
     */
    template<class SCVF>
    static constexpr auto area(const SCVF& scvf)
    {
        static_assert(int(SCVF::Traits::Geometry::mydimension) == int(SCVF::Traits::Geometry::coorddimension-1), "Area element to be called with a codim-1-entity!");
        static_assert(SCVF::Traits::Geometry::coorddimension <= 2, "Axis rotation only makes sense for geometries up to 2D!");
        static_assert(radialAxis < int(SCVF::Traits::Geometry::coorddimension), "Illegal radial axis!");

        // Guldinus theorem
        return scvf.area()*2.0*M_PI*scvf.center()[radialAxis];
    }

    /*!
     * \brief Transformed sub-control-volume volume
     * \note Mid-point rule integrals are only exact for constants
     */
    template<class SCV>
    static constexpr auto volume(const SCV& scv)
    {
        static_assert(int(SCV::Traits::Geometry::mydimension) == int(SCV::Traits::Geometry::coorddimension), "Volume element to be called with a codim-0-entity!");
        static_assert(SCV::Traits::Geometry::coorddimension <= 2, "Axis rotation only makes sense for geometries up to 2D!");
        static_assert(radialAxis < int(SCV::Traits::Geometry::coorddimension), "Illegal radial axis!");

        // Guldinus theorem
        return scv.volume()*2.0*M_PI*scv.center()[radialAxis];
    }

    /*!
     * \brief Integration element for quadrature rules on the reference element
     */
    template<class Geometry>
    static constexpr auto integrationElement(const Geometry& geo, const typename Geometry::LocalCoordinate& x)
    {
        static_assert(Geometry::coorddimension <= 2, "Axis rotation only makes sense for geometries up to 2D!");
        static_assert(radialAxis < int(Geometry::coorddimension), "Illegal radial axis!");

        // Multiply with the polar extrusion factor (2*pi) and the determinant of the transformation Jacobian (radius)
        return geo.integrationElement(x)*2.0*M_PI*geo.global(x)[radialAxis];
    }
};

/*!
 * \ingroup Discretization
 * \brief Rotation symmetric extrusion policy for spherical rotation
 */
struct SphericalExtrusion
{
    /*!
     * \brief Transformed sub-control-volume face area
     * \note Mid-point rule integrals are only exact for constants
     */
    template<class SCVF>
    static constexpr auto area(const SCVF& scvf)
    {
        static_assert(int(SCVF::Traits::Geometry::mydimension) == int(SCVF::Traits::Geometry::coorddimension-1), "Area element to be called with a codim-1-entity!");
        static_assert(SCVF::Traits::Geometry::coorddimension == 1, "Spherical rotation only makes sense for 1D geometries!");

        // sphere surface area
        const auto radius = scvf.center()[0];
        return 4.0*M_PI*radius*radius;
    }

    /*!
     * \brief Transformed sub-control-volume volume
     * \note Mid-point rule integrals are only exact for constants
     */
    template<class SCV>
    static constexpr auto volume(const SCV& scv)
    {
        static_assert(int(SCV::Traits::Geometry::mydimension) == int(SCV::Traits::Geometry::coorddimension), "Volume element to be called with a codim-0-entity!");
        static_assert(SCV::Traits::Geometry::coorddimension == 1, "Spherical rotation only makes sense for 1D geometries!");

        // subtract two balls
        const auto radius0 = scv.corner(0)[0];
        const auto radius1 = scv.corner(1)[0];
        using std::abs;
        return 4.0/3.0*M_PI*abs(radius1*radius1*radius1 - radius0*radius0*radius0);
    }

    /*!
     * \brief Integration element for quadrature rules on the reference element
     */
    template<class Geometry>
    static constexpr auto integrationElement(const Geometry& geo, const typename Geometry::LocalCoordinate& x)
    {
        static_assert(Geometry::coorddimension == 1, "Spherical rotation only makes sense for 1D geometries!");

        // Multiply with the constant spherical extrusion factor (int_0^2pi int_0^pi sin(phi) dphi dtheta = 4pi)
        // and the remaining (radius-dependent) part of determinant of the transformation Jacobian (radius*radius)
        const auto radius = geo.global(x)[0];
        return geo.integrationElement(x)*4.0*M_PI*radius*radius;
    }
};

/*!
 * \brief Traits extracting the public Extrusion type from T
 * Defaults to NoExtrusion if no such type is found
 */
template<class T>
class Extrusion
{
    template<class G>
    using E = typename G::Extrusion;
public:
    using type = typename Dune::Std::detected_or<NoExtrusion, E, T>::type;
};

/*!
 * \brief Convenience alias for obtaining the extrusion type
 */
template<class T>
using Extrusion_t = typename Extrusion<T>::type;

} // end namespace Dumux

#endif
