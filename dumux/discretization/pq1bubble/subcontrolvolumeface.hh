// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the cvfe scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PQ1BubbleDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = PQ1BubbleMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Class for a sub control volume face in the cvfe method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = PQ1BubbleDefaultScvfGeometryTraits<GV> >
class PQ1BubbleSubControlVolumeFace
: public SubControlVolumeFaceBase<PQ1BubbleSubControlVolumeFace<GV, T>, T>
{
    using ThisType = PQ1BubbleSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    PQ1BubbleSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    PQ1BubbleSubControlVolumeFace(const GlobalPosition& center,
                                  const Scalar area,
                                  const GlobalPosition& normal,
                                  const std::array<LocalIndexType, 2>& scvIndices,
                                  const LocalIndexType localScvfIdx,
                                  bool overlapping = false)
    : center_(center)
    , unitOuterNormal_(normal)
    , area_(area)
    , localScvfIdx_(localScvfIdx)
    , scvIndices_(scvIndices)
    , boundary_(false)
    , overlapping_(overlapping)
    , boundaryFlag_{}
    { }

    //! Constructor for boundary scvfs
    PQ1BubbleSubControlVolumeFace(const GlobalPosition& center,
                                  const Scalar area,
                                  const GlobalPosition& normal,
                                  const std::array<LocalIndexType, 2>& scvIndices,
                                  const LocalIndexType localScvfIdx,
                                  const BoundaryFlag& bFlag,
                                  bool overlapping = false)
    : center_(center)
    , unitOuterNormal_(normal)
    , area_(area)
    , localScvfIdx_(localScvfIdx)
    , scvIndices_(scvIndices)
    , boundary_(true)
    , overlapping_(overlapping)
    , boundaryFlag_(bFlag)
    {}

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return center_; }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns true if the sub control volume face is overlapping with another scv
    bool isOverlapping() const
    { return overlapping_; }

    bool boundary() const
    { return boundary_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    { return unitOuterNormal_; }

    //! Index of the inside sub control volume
    GridIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! index of the outside sub control volume
    GridIndexType outsideScvIdx() const
    { return scvIndices_[1]; }

    //! The number of scvs on the outside of this face
    std::size_t numOutsideScvs() const
    { return static_cast<std::size_t>(!boundary()); }

    //! The local index of this sub control volume face
    LocalIndexType index() const
    { return localScvfIdx_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

private:
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    LocalIndexType localScvfIdx_;
    std::array<LocalIndexType, 2> scvIndices_;
    bool boundary_;
    bool overlapping_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux

#endif
