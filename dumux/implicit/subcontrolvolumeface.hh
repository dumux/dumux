// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for a sub control volume face, i.e a part of the boundary
 *        of a sub control volume we computing a flux on.
 */
template<class Geometry, typename IndexType>
class SubControlVolumeFace
{
    using Scalar = typename Geometry::ctype;
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    SubControlVolumeFace(const Geometry& geometry,
                         const GlobalPosition& integrationPoint,
                         const GlobalPosition& unitOuterNormal,
                         IndexType scvfIndex,
                         const std::vector<IndexType>& scvIndices,
                         bool boundary = false)
    : geometry_(geometry),
      integrationPoint_(integrationPoint),
      unitOuterNormal_(unitOuterNormal),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      boundary_(boundary) {}

    //! The center of the sub control volume face
    GlobalPosition center() const
    {
        return geometry_.center();
    }

    //! The integration point for flux evaluations
    GlobalPosition integrationPoint() const
    {
        // Return center for now
        return integrationPoint_;
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return geometry_.volume();
    }

    //! The geometry of the sub control volume face
    const Geometry& geometry() const
    {
        return geometry_;
    }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    GlobalPosition unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume for spatial param evaluation
    IndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is false
    IndexType outsideScvIdx() const
    {
        return scvIndices_[1];
    }

    //! The global index of this sub control volume face
    IndexType index() const
    {
        return scvfIndex_;
    }

private:
    Geometry geometry_;
    GlobalPosition integrationPoint_;
    GlobalPosition unitOuterNormal_;
    IndexType scvfIndex_;
    std::vector<IndexType> scvIndices_;
    bool boundary_;
};

} // end namespace

#endif
