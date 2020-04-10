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
 * \ingroup Common
 * \brief Boundary flag to store e.g. in sub control volume faces
 */
#ifndef DUMUX_BOUNDARY_FLAG_HH
#define DUMUX_BOUNDARY_FLAG_HH

#include <cstddef>
#include <limits>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Class for accessing boundary flags
 * \note this works for all grid managers with gmsh meshes.
 */
class BoundarySegmentIndexFlag
{
public:
    BoundarySegmentIndexFlag()
    : flag_(std::numeric_limits<std::size_t>::max()) {}

    template<class Intersection>
    BoundarySegmentIndexFlag(const Intersection& i)
    : flag_(std::numeric_limits<std::size_t>::max())
    {
        if (i.boundary())
            flag_ = i.boundarySegmentIndex();
    }

    using value_type = std::size_t;

    value_type get() const { return flag_; }

private:
    value_type flag_;
};

/*!
 * \file
 * \ingroup Common
 * \brief Boundary flag to store e.g. in sub control volume faces
 * \note Can be specialized for each grid manager (in the gridmanager headers)
 * \tparam Grid the type of the grid
 */
template<class Grid>
class BoundaryFlag : public BoundarySegmentIndexFlag
{ using BoundarySegmentIndexFlag::BoundarySegmentIndexFlag; };

}  // end namespace Dumux

#endif
