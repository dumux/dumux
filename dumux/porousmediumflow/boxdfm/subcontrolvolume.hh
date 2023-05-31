// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \brief the sub control volume for the box discrete fracture scheme
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUME_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUME_HH

#include <dune/common/reservedvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/porousmediumflow/boxdfm/geometryhelper.hh>
#include <dumux/common/math.hh>
#include <dumux/geometry/volume.hh>

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
    using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
    using LocalIndexType = unsigned int;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = BoxDfmMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim, dimWorld>::Type;
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
    , volume_(Dumux::convexPolytopeVolume<T::dim>(
        Dune::GeometryTypes::cube(T::dim),
        [&](unsigned int i){ return corners_[i]; })
    )
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
        auto corners = geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), indexInIntersection);
        const auto numCorners = corners.size();
        corners_.resize(numCorners);
        for (unsigned int i = 0; i < numCorners; ++i)
            corners_[i] = corners[i];

        // compute volume and scv center
        volume_ = Dumux::convexPolytopeVolume<T::dim-1>(
            Dune::GeometryTypes::cube(T::dim-1),
            [&](unsigned int i){ return corners_[i]; }
        );
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
    [[deprecated("Will be removed after 3.7. Use fvGeometry.geometry(scv).")]]
    Geometry geometry() const
    {
        if (isFractureScv_)
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extract the corners from this scv "
                                                    "and extrude them by the corresponding aperture. ");

        return Geometry(Dune::GeometryTypes::cube(dim), corners_);
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
