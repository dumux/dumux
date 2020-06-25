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
 * \brief Wrapper to make a sub control volume face rotation symmetric
 */
#ifndef DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_SUBCONTROLVOLUMEFACE_HH

#warning "This header is deprecated and will be removed after release 3.3"
#include <cmath>
#include <dumux/discretization/rotationpolicy.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume face rotation symmetric
 * \tparam SubControlVolumeFace The wrapped scvf type
 * \tparam rotationPolicy the rotation policy (see enum RotationPolicy)
 */
template<class SubControlVolumeFace, RotationPolicy rotationPolicy>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolumeFace;

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume face rotation symmetric
 * \tparam SubControlVolumeFace The wrapped scvf type
 * \note Specialization for the 'disc' policy (1d grid --> 2d disc)
 */
template<class SubControlVolumeFace>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolumeFace<SubControlVolumeFace, RotationPolicy::disc>
: public SubControlVolumeFace
{
    using Scalar = typename SubControlVolumeFace::Traits::Scalar;
    static_assert(SubControlVolumeFace::Traits::Geometry::mydimension == 0, "Rotation symmetric scvf with disc policy only works with 1d grids!");
    static_assert(SubControlVolumeFace::Traits::Geometry::coorddimension == 1, "Rotation symmetric scvf with disc policy only works with 1d grids!");
public:
    //! Parent type constructor
    using SubControlVolumeFace::SubControlVolumeFace;

    //! The area of the sub control volume face (circumference of circle)
    Scalar area() const
    {
        const auto radius = this->corner(0)[0];
        return 2.0*M_PI*radius;
    }
};

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume face rotation symmetric
 * \tparam SubControlVolumeFace The wrapped scvf type
 * \note Specialization for the 'ball' policy (1d grid --> 3d ball)
 */
template<class SubControlVolumeFace>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolumeFace<SubControlVolumeFace, RotationPolicy::ball>
: public SubControlVolumeFace
{
    using Scalar = typename SubControlVolumeFace::Traits::Scalar;
    static_assert(SubControlVolumeFace::Traits::Geometry::mydimension == 0, "Rotation symmetric scvf with ball policy only works with 1d grids!");
    static_assert(SubControlVolumeFace::Traits::Geometry::coorddimension == 1, "Rotation symmetric scvf with ball policy only works with 1d grids!");
public:
    //! Parent type constructor
    using SubControlVolumeFace::SubControlVolumeFace;

    //! The area of the sub control volume face (surface of sphere)
    Scalar area() const
    {
        const auto radius = this->corner(0)[0];
        return 4.0*M_PI*radius*radius;
    }
};

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume face rotation symmetric
 * \tparam SubControlVolumeFace The wrapped scvf type
 * \note Specialization for the 'toroid' policy (2d grid --> 3d toroid)
 * \note We rotate about the second axis
 */
template<class SubControlVolumeFace>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolumeFace<SubControlVolumeFace, RotationPolicy::toroid>
: public SubControlVolumeFace
{
    using Scalar = typename SubControlVolumeFace::Traits::Scalar;
    static_assert(SubControlVolumeFace::Traits::Geometry::mydimension == 1, "Rotation symmetric scvf with toroid policy only works with 2d grids!");
    static_assert(SubControlVolumeFace::Traits::Geometry::coorddimension == 2, "Rotation symmetric scvf with toroid policy only works with 2d grids!");
public:
    //! Parent type constructor
    using SubControlVolumeFace::SubControlVolumeFace;

    //! The area of the sub control volume face (Guldinus theorem)
    Scalar area() const
    {
        const auto radius = this->center()[0];
        return SubControlVolumeFace::area()*2.0*M_PI*radius;
    }
};

} // end namespace Dumux

#endif
