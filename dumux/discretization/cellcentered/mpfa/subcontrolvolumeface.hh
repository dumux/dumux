// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief The sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <vector>
#include <array>

#include <dune/common/reservedvector.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/boundaryflag.hh>

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
            static const unsigned int topologyId = Dune::GeometryTypes::cube(dim).id();
        };
    };

    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar> >;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
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
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    // Information on the intersection from which this scvf was constructed
    struct FacetInfo
    {
        GridIndexType elementIndex;
        int facetIndex;
        int facetCornerIndex;
    };

    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    /*!
     * \brief Constructor
     *
     * \param helper The helper class for mpfa schemes
     * \param corners The corners of the scv face
     * \param intersection The intersection
     * \param facetInfo Information on the facet from which this scvf is constructed
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param vIdxLocal The element-local vertex index the scvf is connected to
     * \param scvfIndex The global index of this scv face
     * \param insideScvIdx The inside scv index connected to this face
     * \param outsideScvIndices The outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    //! Construction with given intersection
    template<class MpfaHelper, class Intersection>
    CCMpfaSubControlVolumeFace(const MpfaHelper& helper,
                               CornerStorage&& corners,
                               const Intersection& intersection,
                               FacetInfo facetInfo,
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
    , center_(0.0)
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , boundaryFlag_{intersection}
    , facetInfo_{std::move(facetInfo)}
    {
          // compute the center of the scvf
          for (const auto& corner : corners)
              center_ += corner;
          center_ /= corners.size();

          // use helper class to obtain area & integration point
          ipGlobal_ = helper.getScvfIntegrationPoint(corners, q);
          area_ = helper.computeScvfArea(corners);
    }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns true if the sub control volume face is on the domain boundary
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

    //! The number of scvs on the outside of this scv face
    std::size_t numOutsideScvs() const
    { return outsideScvIndices_.size(); }

    //! Index of the i-th outside sub control volume or boundary scv index.
    // Results in undefined behaviour if i >= numOutsideScvs()
    GridIndexType outsideScvIdx(int i = 0) const
    { return outsideScvIndices_[i]; }

    //! returns the outside scv indices (can be more than one index for dim < dimWorld)
    const OutsideGridIndexStorage& outsideScvIndices() const
    { return outsideScvIndices_; }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! returns the unit outer normal vector (assumes non-curved geometries)
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

    //! Return information on the facet from which this scvf was constructed
    const FacetInfo& facetInfo() const
    { return facetInfo_; }

private:
    bool boundary_;
    GridIndexType vertexIndex_;
    GridIndexType scvfIndex_;
    GridIndexType insideScvIdx_;
    OutsideGridIndexStorage outsideScvIndices_;
    unsigned int vIdxInElement_;

    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    BoundaryFlag boundaryFlag_;
    FacetInfo facetInfo_;
};

} // end namespace Dumux

#endif
