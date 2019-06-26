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
 * \ingroup CCMpfaDiscretization
 * \brief The sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH

#include <vector>
#include <array>

#include <dune/common/reservedvector.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the cell-centered finite volume scheme using MPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct CCMpfaDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using Scalar = typename Grid::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using OutsideGridIndexStorage = typename std::conditional_t< (dim<dimWorld),
                                                                 std::vector<GridIndexType>,
                                                                 Dune::ReservedVector<GridIndexType, 1> >;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int dim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< dim >::type::id;
        };
    };

    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar> >;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for a sub control volume face in mpfa methods, i.e a part of the boundary
 *        of a control volume we compute fluxes on.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = CCMpfaDefaultScvfGeometryTraits<GV> >
class CCMpfaSubControlVolumeFace
{
    using GridIndexType = typename T::GridIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using OutsideGridIndexStorage = typename T::OutsideGridIndexStorage;
    using Geometry = typename T::Geometry;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    /*!
     * \brief Constructor
     *
     * \param helper The helper class for mpfa schemes
     * \param corners The corners of the scv face
     * \param unitOuterNormal The unit outer normal vector of the scvf
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param vIdxLocal The element-local vertex index the scvf is connected to
     * \param scvfIndex The global index of this scv face
     * \param insideScvIdx The inside scv index connected to this face
     * \param outsideScvIndices The outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    template<class MpfaHelper>
    CCMpfaSubControlVolumeFace(const MpfaHelper& helper,
                               CornerStorage&& corners,
                               GlobalPosition&& unitOuterNormal,
                               GridIndexType vIdxGlobal,
                               unsigned int vIdxLocal,
                               GridIndexType scvfIndex,
                               GridIndexType insideScvIdx,
                               const OutsideGridIndexStorage& outsideScvIndices,
                               Scalar q,
                               bool boundary)
    : boundary_(boundary)
    , vertexIndex_(vIdxGlobal)
    , scvfIndex_(scvfIndex)
    , insideScvIdx_(insideScvIdx)
    , outsideScvIndices_(outsideScvIndices)
    , vIdxInElement_(vIdxLocal)
    , corners_(std::move(corners))
    , center_(0.0)
    , unitOuterNormal_(std::move(unitOuterNormal))
    {
          // compute the center of the scvf
          for (const auto& corner : corners_)
              center_ += corner;
          center_ /= corners_.size();

          // use helper class to obtain area & integration point
          ipGlobal_ = helper.getScvfIntegrationPoint(corners_, q);
          area_ = helper.computeScvfArea(corners_);
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

    //! Returns the element-local vertex index the scvf is connected to
    unsigned int vertexIndexInElement() const
    { return vIdxInElement_; }

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
    const OutsideGridIndexStorage& outsideScvIndices() const
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
    OutsideGridIndexStorage outsideScvIndices_;
    unsigned int vIdxInElement_;

    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
};

} // end namespace Dumux

#endif
