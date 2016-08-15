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
 * \brief Class for an MPFA-O sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 */
template<class G, typename I>
class CCMpfaSubControlVolumeFace : public SubControlVolumeFaceBase<CCMpfaSubControlVolumeFace<G, I>, G, I>
{
    using ParentType = SubControlVolumeFaceBase<CCMpfaSubControlVolumeFace<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;

public:
    template<class MpfaGeometryHelper>
    CCMpfaSubControlVolumeFace(const MpfaGeometryHelper& geomHelper,
                               std::vector<GlobalPosition>&& corners,
                               GlobalPosition&& unitOuterNormal,
                               IndexType vertexIndex,
                               IndexType scvfIndex,
                               std::array<IndexType, 2>&& scvIndices,
                               Scalar q,
                               bool boundary)
    : ParentType(),
      boundary_(boundary),
      vertexIndex_(vertexIndex),
      scvfIndex_(scvfIndex),
      scvIndices_(std::move(scvIndices)),
      corners_(std::move(corners)),
      center_(0.0),
      unitOuterNormal_(std::move(unitOuterNormal))
      {
            for (const auto& corner : corners_)
                center_ += corner;
            center_ /= corners_.size();
            ipGlobal_ = geomHelper.getScvfIntegrationPoint(corners_, q);
            area_ = geomHelper.getScvfArea(corners_);
      }

    //! The center of the sub control volume face
    GlobalPosition center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    GlobalPosition ipGlobal() const
    { return ipGlobal_; }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    { return boundary_; }

    GlobalPosition unitOuterNormal() const
    { return unitOuterNormal_; }

    //! index of the inside sub control volume for spatial param evaluation
    IndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! index of the outside sub control volume for spatial param evaluation
    IndexType outsideScvIdx() const
    { return scvIndices_[1]; }

    //! The global index of this sub control volume face
    IndexType index() const
    { return scvfIndex_; }

    GlobalPosition corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    const Geometry geometry() const
    { return Geometry(Dune::GeometryType(Dune::GeometryType::cube, dim), corners_); }

    //! Returns the global position of the vertex the scvf is connected to
    GlobalPosition vertexCorner() const
    { return corner(corners_.size()-1); }

    //! Returns the global position of the center of the element facet this scvf is embedded in
    GlobalPosition facetCorner() const
    { return corner(0); }

    //! Returns the index of the vertex the scvf is connected to
    IndexType vertexIndex() const
    { return vertexIndex_; }

private:
    bool boundary_;
    IndexType vertexIndex_;
    IndexType scvfIndex_;
    std::array<IndexType, 2> scvIndices_;

    std::vector<GlobalPosition> corners_;
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
};

} // end namespace

#endif
