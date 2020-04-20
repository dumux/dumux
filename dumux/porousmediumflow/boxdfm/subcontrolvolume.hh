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
 * \ingroup BoxDFMModel
 * \brief the sub control volume for the box discrete fracture scheme
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUME_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUME_HH

#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief Default traits class to be used for the sub-control volumes
 *        for the box discrete fracture scheme
 *
 * \tparam GV the type of the grid view
 *
 * \note We define new traits for the box-dfm sub-control volume face
 *       as we use a different type of container for storing the scvf corners!
 */
template<class GridView>
struct BoxDfmDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    template <class ct>
    struct ScvMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim) corners (1<<(dim))
        // However, on fracture scvs the number might be smaller (use ReservedVector)
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<(dim)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
    using LocalIndexType = unsigned int;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, ScvMLGTraits<Scalar>>;
    using CornerStorage = typename ScvMLGTraits<Scalar>::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
};

/*!
 * \ingroup BoxDFMModel
 * \brief the sub control volume for the box discrete fracture scheme
 *
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDfmDefaultScvGeometryTraits<GV> >
class BoxDfmSubControlVolume
: public SubControlVolumeBase<BoxDfmSubControlVolume<GV, T>, T>
{
    using ThisType = BoxDfmSubControlVolume<GV, T>;
    using ParentType = SubControlVolumeBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using GlobalPosition = typename T::GlobalPosition;
    using CornerStorage = typename T::CornerStorage;
    enum { dim = Geometry::mydimension };

    static_assert(dim == 2 || dim == 3, "Box-Dfm sub-control volume only implemented in 2d or 3d");

public:
    //! State the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxDfmSubControlVolume() = default;

    // the constructor for standard scvs
    template<class GeometryHelper>
    BoxDfmSubControlVolume(const GeometryHelper& geometryHelper,
                           LocalIndexType scvIdx,
                           GridIndexType elementIndex,
                           GridIndexType dofIndex)
    : isFractureScv_(false)
    , corners_(geometryHelper.getScvCorners(scvIdx))
    , center_(0.0)
    , volume_(geometryHelper.scvVolume(corners_))
    , elementIndex_(elementIndex)
    , vIdxLocal_(scvIdx)
    , elemLocalScvIdx_(scvIdx)
    , dofIndex_(dofIndex)
    , facetIdx_(0)
    {
        // compute center point
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    /*!
     * \brief Constructor for fracture scvs
     *
     * The corner computation is the same as for boundary scvfs.
     * Also, the scvf area of a boundary scvf is equal to the scv
     * volume (unscaled by the aperture) Thus, we reuse functionality here.
     * In order to get the right dimensions later, one must provide appropriate
     * extrusion factors in the problem corresponding to the fracture aperture.     *
     */
    template<class GeometryHelper, class Intersection>
    BoxDfmSubControlVolume(const GeometryHelper& geometryHelper,
                           const Intersection& intersection,
                           const typename Intersection::Geometry& isGeometry,
                           LocalIndexType indexInIntersection,
                           LocalIndexType vIdxLocal,
                           LocalIndexType elemLocalScvIdx,
                           LocalIndexType elemLocalFacetIdx,
                           GridIndexType elementIndex,
                           GridIndexType dofIndex)
    : isFractureScv_(true)
    , corners_()
    , center_(0.0)
    , volume_(0.0)
    , elementIndex_(elementIndex)
    , vIdxLocal_(vIdxLocal)
    , elemLocalScvIdx_(elemLocalScvIdx)
    , dofIndex_(dofIndex)
    , facetIdx_(elemLocalFacetIdx)
    {
        // copy corners
        auto corners = geometryHelper.getBoundaryScvfCorners(intersection, isGeometry, indexInIntersection);
        const auto numCorners = corners.size();
        corners_.resize(numCorners);
        for (unsigned int i = 0; i < numCorners; ++i)
            corners_[i] = corners[i];

        // compute volume and scv center
        volume_ = geometryHelper.scvfArea(corners);
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The volume of the sub control volume
    Scalar volume() const
    { return volume_; }

    //! The element-local vertex index this scv is connected to
    LocalIndexType localDofIndex() const
    { return vIdxLocal_; }

    //! The element-local index of this scv
    LocalIndexType indexInElement() const
    { return elemLocalScvIdx_; }

    //! The element-local facet index for which a fracture scv was created
    LocalIndexType facetIndexInElement() const
    { assert(isFractureScv_); return facetIdx_; }

    //! The index of the dof this scv is embedded in
    GridIndexType dofIndex() const
    { return dofIndex_; }

    // The position of the dof this scv is embedded in (list is defined such that first entry is vertex itself)
    const GlobalPosition& dofPosition() const
    { return corners_[0]; }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    { return elementIndex_; }

    //! Return true if this scv is part of the fracture domain
    bool isOnFracture() const
    { return isFractureScv_; }

    //! The geometry of the sub control volume
    // e.g. for integration
    Geometry geometry() const
    {
        if (isFractureScv_)
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extract the corners from this scv "
                                                    "and extrude them by the corresponding aperture. ");

        return Geometry(Dune::GeometryTypes::cube(dim), corners_);
    }

    //! Return the corner for the given local index
    const GlobalPosition& corner(LocalIndexType localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

private:
    bool isFractureScv_;
    CornerStorage corners_;
    GlobalPosition center_;
    Scalar volume_;
    GridIndexType elementIndex_;
    LocalIndexType vIdxLocal_;
    LocalIndexType elemLocalScvIdx_;
    GridIndexType dofIndex_;

    // for fracture scvs only!
    LocalIndexType facetIdx_;
};

} // end namespace

#endif
