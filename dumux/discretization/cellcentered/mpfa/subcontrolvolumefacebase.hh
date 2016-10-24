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
 * \brief Class for an mpfa-o sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACEBASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACEBASE_HH

#include <utility>
#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Base class for a sub-control volume face in mpfa methods. All mpfa method-specific implementations should inherit from this class
 */
template<class G, typename I>
class CCMpfaSubControlVolumeFaceBase : public SubControlVolumeFaceBase<CCMpfaSubControlVolumeFaceBase<G, I>, G, I>
{
    using ParentType = SubControlVolumeFaceBase<CCMpfaSubControlVolumeFaceBase<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    /*!
     * \brief Constructor
     *
     * We do not use the localIndex here. Its meaning can vary depending on the
     * implementation (i.e. mpfa method) and is handled by the implementation itself.
     *
     * \param geomHelper The mpfa geometry helper
     * \param corners The corners of the scv face
     * \param unitOuterNormal The unit outer normal vector of the scvf
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param localIndex Some element local index (the local vertex index in mpfao-fps)
     * \param scvfIndex The global index of this scv face
     * \param scvIndices The inside and outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    template<class MpfaGeometryHelper>
    CCMpfaSubControlVolumeFaceBase(const MpfaGeometryHelper& geomHelper,
                                   std::vector<GlobalPosition>&& corners,
                                   GlobalPosition&& unitOuterNormal,
                                   IndexType vertexIndex,
                                   unsigned int localIndex,
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

    //! returns bolean if the sub control volume face is on the domain boundary
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

    //! Returns the number of corners
    std::size_t corners() const
    { return corners_.size(); }

    //! Returns the corner for a given local index
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
    { return corners_.back(); }

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
