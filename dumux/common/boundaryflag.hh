// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Boundary flag to store e.g. in sub control volume faces
 */
#ifndef DUMUX_BOUNDARY_FLAG_HH
#define DUMUX_BOUNDARY_FLAG_HH

#include <cstddef>
#include <limits>

namespace Dumux {

/*!
 * \ingroup Core
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
 * \ingroup Core
 * \brief Boundary flag to store e.g. in sub control volume faces
 * \note Can be specialized for each grid manager (in the gridmanager headers)
 * \tparam Grid the type of the grid
 */
template<class Grid>
class BoundaryFlag : public BoundarySegmentIndexFlag
{ using BoundarySegmentIndexFlag::BoundarySegmentIndexFlag; };

}  // end namespace Dumux

#endif
