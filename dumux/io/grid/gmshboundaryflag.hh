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
 * \ingroup InputOutput
 * \brief Boundary flag implementation for Gmsh meshes.
 */

#ifndef DUMUX_GMSH_BOUNDARY_FLAG_HH
#define DUMUX_GMSH_BOUNDARY_FLAG_HH


namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Class for accessing boundary flags for Gmsh meshes.
 */
template<class Grid>
class GmshBoundaryFlag
{
public:
    GmshBoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    GmshBoundaryFlag(const Intersection& i) : flag_(-1)
    {
        if (i.boundary())
            flag_ = i.boundarySegmentIndex();
    }

    using value_type = std::size_t;

    value_type get() const { return flag_; }

private:
    value_type flag_;
};

} // end namespace Dumux

#endif
