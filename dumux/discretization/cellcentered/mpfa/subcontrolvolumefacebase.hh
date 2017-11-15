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
#include <dune/geometry/type.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Base class for a sub-control volume face in mpfa methods.
 *        All mpfa method-specific implementations should inherit from this class
 */
template<class ScvfGeometryTraits>
class CCMpfaSubControlVolumeFaceBase : public SubControlVolumeFaceBase<CCMpfaSubControlVolumeFaceBase<ScvfGeometryTraits>, ScvfGeometryTraits>
{
    using ParentType = SubControlVolumeFaceBase<CCMpfaSubControlVolumeFaceBase<ScvfGeometryTraits>, ScvfGeometryTraits>;
    using GridIndexType = typename ScvfGeometryTraits::GridIndexType;
    using Scalar = typename ScvfGeometryTraits::Scalar;
    using GlobalPosition = typename ScvfGeometryTraits::GlobalPosition;
    using CornerStorage = typename ScvfGeometryTraits::CornerStorage;
    using Geometry = typename ScvfGeometryTraits::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = ScvfGeometryTraits;

    /*!
     * \brief Constructor
     *
     * \param geomHelper The mpfa geometry helper
     * \param corners The corners of the scv face
     * \param unitOuterNormal The unit outer normal vector of the scvf
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param scvfIndex The global index of this scv face
     * \param insideScvIdx The inside scv index connected to this face
     * \param outsideScvIndices The outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    template<class MpfaHelper>
    CCMpfaSubControlVolumeFaceBase(const MpfaHelper& helper,
                                   CornerStorage&& corners,
                                   GlobalPosition&& unitOuterNormal,
                                   GridIndexType vIdxGlobal,
                                   GridIndexType scvfIndex,
                                   GridIndexType insideScvIdx,
                                   const std::vector<GridIndexType>& outsideScvIndices,
                                   Scalar q,
                                   bool boundary)
    : ParentType(),
      boundary_(boundary),
      vertexIndex_(vIdxGlobal),
      scvfIndex_(scvfIndex),
      insideScvIdx_(insideScvIdx),
      outsideScvIndices_(outsideScvIndices),
      corners_(std::move(corners)),
      center_(0.0),
      unitOuterNormal_(std::move(unitOuterNormal))
      {
            // compute the center of the scvf
            for (const auto& corner : corners_)
                center_ += corner;
            center_ /= corners_.size();

            // use helper class to obtain area & integration point
            ipGlobal_ = helper.getScvfIntegrationPoint(corners_, q);
            area_ = helper.getScvfArea(corners_);
      }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns bolean if the sub control volume face is on the domain boundary
    bool boundary() const
    { return boundary_; }

    //! The global index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

    //! Returns the index of the vertex the scvf is connected to
    GridIndexType vertexIndex() const
    { return vertexIndex_; }

    //! index of the inside sub control volume
    GridIndexType insideScvIdx() const
    { return insideScvIdx_; }

    //! The number of outside scvs connection via this scv face
    std::size_t numOutsideScvs() const
    { return outsideScvIndices_.size(); }

    //! index of the outside sub control volume or boundary scv index
    //! returns undefined behaviour if index exceeds numOutsideScvs
    GridIndexType outsideScvIdx(int i = 0) const
    { return outsideScvIndices_[i]; }

    //! returns the outside scv indices (can be more than one index for dim < dimWorld)
    const std::vector<GridIndexType>& outsideScvIndices() const
    { return outsideScvIndices_; }

    //! Returns the number of corners
    std::size_t corners() const
    { return corners_.size(); }

    //! Returns the corner for a given local index
    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! Returns the global position of the vertex the scvf is connected to
    const GlobalPosition& vertexCorner() const
    { return corners_.back(); }

    //! Returns the global position of the center of the element facet this scvf is embedded in
    const GlobalPosition& facetCorner() const
    { return corner(0); }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! returns the unit outer normal vector (assumes non-curved geometries)
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    { return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_); }

private:
    bool boundary_;
    GridIndexType vertexIndex_;
    GridIndexType scvfIndex_;
    GridIndexType insideScvIdx_;
    std::vector<GridIndexType> outsideScvIndices_;

    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
};

} // end namespace

#endif
