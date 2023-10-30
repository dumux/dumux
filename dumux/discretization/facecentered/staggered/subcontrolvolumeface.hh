// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Default traits class to be used for the sub-control volume face
 *        for the face-centered staggered finite volume scheme
 * \tparam GridView the type of the grid view
 */
template<class GridView>
struct FaceCenteredDefaultScvfGeometryTraits
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::Grid::dimension;
    static constexpr int dimWorld = GridView::Grid::dimensionworld;
    using CornerStorage = std::array<GlobalPosition, (1<<(dim-1))>;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim-1, dimWorld>;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Face centered staggered sub control volume face
 */
template<class GridView, class T = FaceCenteredDefaultScvfGeometryTraits<GridView>>
class FaceCenteredStaggeredSubControlVolumeFace
{
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using Scalar = typename T::Scalar;
    using Element = typename T::Element;
    using CornerStorage = typename T::CornerStorage;

    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    using ElementGeometry = typename Element::Geometry;
    using IntersectionGeometry = typename GridView::Intersection::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    using GlobalPosition = typename T::GlobalPosition;
    enum class FaceType : SmallLocalIndexType {frontal, lateral};
    enum class BoundaryType : SmallLocalIndexType {interior, physicalBoundary, processorBoundary};

    FaceCenteredStaggeredSubControlVolumeFace() = default;

    //! The constructor for frontal faces
    FaceCenteredStaggeredSubControlVolumeFace(const ElementGeometry& elementGeometry,
                                              const IntersectionGeometry& intersectionGeometry,
                                              const std::array<GridIndexType, 2> globalScvIndices,
                                              const SmallLocalIndexType localScvfIdx,
                                              const GridIndexType globalScvfIdx,
                                              const GlobalPosition& unitOuterNormal,
                                              const FaceType faceType,
                                              const BoundaryType boundaryType)
    : globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , globalScvfIdx_(globalScvfIdx)
    , area_(intersectionGeometry.volume())
    , normalAxis_(Dumux::normalAxis(unitOuterNormal))
    , outerNormalSign_(sign(unitOuterNormal[normalAxis_]))
    , faceType_(faceType)
    , boundaryType_(boundaryType)
    {
        assert(faceType == FaceType::frontal);
        center_ = boundary() ? intersectionGeometry.center() : elementGeometry.center();
        ipGlobal_ = center_;

        if (!boundary())
            outerNormalSign_ *= -1.0;
    }

    //! The constructor for lateral faces
    template<class LateralFacetGeometry>
    FaceCenteredStaggeredSubControlVolumeFace(const ElementGeometry& elementGeometry,
                                              const IntersectionGeometry& intersectionGeometry,
                                              const LateralFacetGeometry& lateralFacetGeometry,
                                              const std::array<GridIndexType, 2> globalScvIndices,
                                              const SmallLocalIndexType localScvfIdx,
                                              const GridIndexType globalScvfIdx,
                                              const GlobalPosition& unitOuterNormal,
                                              const FaceType faceType,
                                              const BoundaryType boundaryType)
    : globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , globalScvfIdx_(globalScvfIdx)
    , area_(0.5*lateralFacetGeometry.volume())
    , normalAxis_(Dumux::normalAxis(unitOuterNormal))
    , outerNormalSign_(sign(unitOuterNormal[normalAxis_]))
    , faceType_(faceType)
    , boundaryType_(boundaryType)
    {
        assert(faceType == FaceType::lateral);
        const auto shift = intersectionGeometry.center() - elementGeometry.center();
        ipGlobal_ = lateralFacetGeometry.center() + shift;
        center_ = 0.5*(lateralFacetGeometry.center() + ipGlobal_);
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point of the sub control volume face
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    {
        GlobalPosition result(0.0);
        result[normalAxis_] = 1.0 * directionSign();
        return result;
    }

    //! Index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    { return globalScvIndices_[0]; }

    //! index of the outside sub control volume for spatial param evaluation
    GridIndexType outsideScvIdx() const
    { return globalScvIndices_[1]; }

    GridIndexType index() const
    { return globalScvfIdx_; }

    SmallLocalIndexType localIndex() const
    { return localScvfIdx_; }

    FaceType faceType() const
    { return faceType_; }

    bool boundary() const
    { return boundaryType_ == BoundaryType::physicalBoundary; }

    bool processorBoundary() const
    { return boundaryType_ == BoundaryType::processorBoundary; }

    bool isFrontal() const
    { return faceType_ == FaceType::frontal; }

    bool isLateral() const
    { return faceType_ == FaceType::lateral; }

    Scalar area() const
    { return area_; }

    SmallLocalIndexType normalAxis() const
    { return normalAxis_; }

    std::int_least8_t directionSign() const
    { return outerNormalSign_; }

private:
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    std::array<GridIndexType, 2> globalScvIndices_;
    SmallLocalIndexType localScvfIdx_;
    GridIndexType globalScvfIdx_;
    Scalar area_;
    SmallLocalIndexType normalAxis_;
    std::int_least8_t outerNormalSign_;
    FaceType faceType_;
    BoundaryType boundaryType_;
};

} // end namespace Dumux

#endif
