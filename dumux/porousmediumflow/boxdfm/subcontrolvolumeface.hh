// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \brief The sub control volume face class for the box discrete fracture model.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/porousmediumflow/boxdfm/geometryhelper.hh>
#include <dumux/geometry/volume.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the box discrete fracture scheme
 *
 * \tparam GV the type of the grid view
 *
 * \note We define new traits for the box-dfm sub-control volume face
 *       as we use a different type of container for storing the scvf corners!
 */
template<class GridView>
struct BoxDfmDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
    using LocalIndexType = unsigned int;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = BoxDfmMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup BoxDFMModel
 * \brief Class for a sub control volume face in the box discrete fracture method, i.e a
 *        part of the boundary of a sub control volume we compute fluxes on.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDfmDefaultScvfGeometryTraits<GV> >
class BoxDfmSubControlVolumeFace
: public SubControlVolumeFaceBase<BoxDfmSubControlVolumeFace<GV, T>, T>
{
    using ThisType = BoxDfmSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using GlobalPosition = typename T::GlobalPosition;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

    static_assert(T::dim == 2 || T::dim == 3, "Box-Dfm sub-control volume face only implemented in 2d or 3d");

public:
    //! State the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxDfmSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class GeometryHelper, class Element>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Element& element,
                               const typename Element::Geometry& elemGeometry,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices)
    : corners_(geometryHelper.getScvfCorners(scvfIndex))
    , center_(0.0)
    , unitOuterNormal_(geometryHelper.normal(corners_, scvIndices))
    , area_(Dumux::convexPolytopeVolume<T::dim-1>(
        Dune::GeometryTypes::cube(T::dim-1),
        [&](unsigned int i){ return corners_[i]; })
    )
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(false)
    , isFractureScvf_(false)
    , boundaryFlag_{}
    , facetIdx_(0)
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for boundary scvfs
    template<class GeometryHelper, class Intersection>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Intersection& intersection,
                               const typename Intersection::Geometry& isGeometry,
                               LocalIndexType indexInIntersection,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices)
    : corners_(geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), indexInIntersection))
    , center_(0.0)
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , area_(Dumux::convexPolytopeVolume<T::dim-1>(
        Dune::GeometryTypes::cube(T::dim-1),
        [&](unsigned int i){ return corners_[i]; })
    )
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(true)
    , isFractureScvf_(false)
    , boundaryFlag_{intersection}
    , facetIdx_(0)
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for inner fracture scvfs
    template<class GeometryHelper, class Intersection>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Intersection& intersection,
                               const typename Intersection::Geometry& isGeometry,
                               LocalIndexType indexInIntersection,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices,
                               bool boundary)
    : corners_(geometryHelper.getFractureScvfCorners(intersection.indexInInside(), indexInIntersection))
    , center_(0.0)
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(boundary)
    , isFractureScvf_(true)
    , boundaryFlag_{intersection}
    , facetIdx_(intersection.indexInInside())
    {
        // The area here is given in meters. In order to
        // get the right dimensions, the user has to provide
        // the appropriate aperture in the problem (via an extrusion factor)
        if (T::dim == 3)
            area_ = (corners_[1]-corners_[0]).two_norm();
        else if (T::dim == 2)
            area_ = 1.0;

        // obtain the unit normal vector
        unitOuterNormal_ = geometryHelper.fractureNormal(corners_, intersection, indexInIntersection);

        // compute the scvf center
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return center_; }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns true if the sub control volume face is on the boundary
    bool boundary() const
    { return boundary_; }

    //! returns the unit normal vector pointing outwards
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! The global index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

    //! Return if this is a fracture scvf
    bool isOnFracture() const
    { return isFractureScvf_; }

    //! The element-local facet index for which a fracture scv was created
    LocalIndexType facetIndexInElement() const
    { assert(isFractureScvf_); return facetIdx_; }

    //! Returns the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

    //! index of the inside sub control volume
    LocalIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! Index of the i-th outside sub control volume or boundary scv index.
    // Results in undefined behaviour if i >= numOutsideScvs()
    LocalIndexType outsideScvIdx(int i = 0) const
    {
        assert(!boundary());
        return scvIndices_[1];
    }

    //! The number of scvs on the outside of this face
    std::size_t numOutsideScvs() const
    {
        return static_cast<std::size_t>(!boundary());
    }

    //! The geometry of the sub control volume face
    [[deprecated("Will be removed after 3.7. Use fvGeometry.geometry(scvf).")]]
    Geometry geometry() const
    {
        if (isFractureScvf_)
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extract the corners from this scv "
                                                    "and extrude them by the corresponding aperture. ");

        return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_);
    }

private:
    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;
    bool boundary_;
    bool isFractureScvf_;
    BoundaryFlag boundaryFlag_;
    LocalIndexType facetIdx_;
};

} // end namespace Dumux

#endif
