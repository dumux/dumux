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
 * \brief Wrapper to make a sub control volume rotation symmetric
 */
#ifndef DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_SUBCONTROLVOLUME_HH

#warning "This header is deprecated and will be removed after release 3.3"
#include <cmath>
#include <dumux/discretization/rotationpolicy.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume rotation symmetric
 * \tparam SubControlVolume The wrapped scv type
 * \tparam rotationPolicy the rotation policy (see enum RotationPolicy)
 */
template<class SubControlVolume, RotationPolicy rotationPolicy>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolume;

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume rotation symmetric
 * \tparam SubControlVolume The wrapped scv type
 * \note Specialization for the 'disc' policy (1d grid --> 2d disc)
 */
template<class SubControlVolume>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolume<SubControlVolume, RotationPolicy::disc>
: public SubControlVolume
{
    using Scalar = typename SubControlVolume::Traits::Scalar;
    static_assert(SubControlVolume::Traits::Geometry::mydimension == 1, "Rotation symmetric scv with disc policy only works with 1d grids!");
    static_assert(SubControlVolume::Traits::Geometry::coorddimension == 1, "Rotation symmetric scv with disc policy only works with 1d grids!");
public:
    //! Parent type constructor
    using SubControlVolume::SubControlVolume;

    //! The volume of the sub control volume (difference between two disks)
    Scalar volume() const
    {
        using std::abs;
        const auto radius0 = this->corner(0)[0];
        const auto radius1 = this->corner(1)[0];
        return M_PI*abs(radius1*radius1 - radius0*radius0);
    }
};

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume rotation symmetric
 * \tparam SubControlVolume The wrapped scv type
 * \note Specialization for the 'ball' policy (1d grid --> 3d ball)
 */
template<class SubControlVolume>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolume<SubControlVolume, RotationPolicy::ball>
: public SubControlVolume
{
    using Scalar = typename SubControlVolume::Traits::Scalar;
    static_assert(SubControlVolume::Traits::Geometry::mydimension == 1, "Rotation symmetric scv with ball policy only works with 1d grids!");
    static_assert(SubControlVolume::Traits::Geometry::coorddimension == 1, "Rotation symmetric scv with ball policy only works with 1d grids!");
public:
    //! Parent type constructor
    using SubControlVolume::SubControlVolume;

    //! The volume of the sub control volume (difference between two balls)
    Scalar volume() const
    {
        using std::abs;
        const auto radius0 = this->corner(0)[0];
        const auto radius1 = this->corner(1)[0];
        return 4.0/3.0*M_PI*abs(radius1*radius1*radius1 - radius0*radius0*radius0);
    }
};

/*!
 * \ingroup Discretization
 * \brief Wrapper to make a sub control volume rotation symmetric
 * \tparam SubControlVolume The wrapped scv type
 * \note Specialization for the 'toroid' policy (2d grid --> 3d toroid)
 * \note We rotate about the second axis
 */
template<class SubControlVolume>
class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricSubControlVolume<SubControlVolume, RotationPolicy::toroid>
: public SubControlVolume
{
    using Scalar = typename SubControlVolume::Traits::Scalar;
    static_assert(SubControlVolume::Traits::Geometry::mydimension == 2, "Rotation symmetric scv with toroid policy only works with 2d grids!");
    static_assert(SubControlVolume::Traits::Geometry::coorddimension == 2, "Rotation symmetric scv with toroid policy only works with 2d grids!");
public:
    //! Parent type constructor
    using SubControlVolume::SubControlVolume;

    //! The volume of the sub control volume (Guldinus theorem)
    Scalar volume() const
    {
        const auto radius = this->center()[0];
        return SubControlVolume::volume()*2.0*M_PI*radius;
    }
};

} // end namespace Dumux

#endif
