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
 * \ingroup PQ2Discretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_PQ2_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/pq2/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the cvfe scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PQ2DefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = PQ2MLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup PQ2Discretization
 * \brief Class for a sub control volume face in the cvfe method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = PQ2DefaultScvfGeometryTraits<GV> >
class PQ2SubControlVolumeFace
: public SubControlVolumeFaceBase<PQ2SubControlVolumeFace<GV, T>, T>
{
    using ThisType = PQ2SubControlVolumeFace<GV, T>;
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
    PQ2SubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    PQ2SubControlVolumeFace(const GlobalPosition& center,
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
    PQ2SubControlVolumeFace(const GlobalPosition& center,
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
