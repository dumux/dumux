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
 * \brief the sub control volume for the cvfe scheme
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_PQ2_SUBCONTROLVOLUME_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/discretization/pq2/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the pq2 scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PQ2DefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = PQ2MLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
};

/*!
 * \ingroup PQ2Discretization
 * \brief the sub control volume for the pq2 scheme
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GridView, class T = PQ2DefaultScvGeometryTraits<GridView>>
class PQ2SubControlVolume
{
    using GlobalPosition = typename T::GlobalPosition;
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    PQ2SubControlVolume() = default;

    PQ2SubControlVolume(const Scalar& volume,
                              const GlobalPosition& dofPosition,
                              const GlobalPosition& center,
                              const LocalIndexType indexInElement,
                              const GridIndexType eIdx,
                              const GridIndexType dofIdx,
                              bool overlapping = false)

    : center_(center)
    , dofPosition_(dofPosition)
    , volume_(volume)
    , indexInElement_(indexInElement)
    , eIdx_(eIdx)
    , dofIdx_(dofIdx)
    , overlapping_(overlapping)
    {}

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The position of the degree of freedom
    const GlobalPosition& dofPosition() const
    { return dofPosition_; }

    Scalar volume() const
    { return volume_; }

    //! returns true if the sub control volume is overlapping with another scv
    bool isOverlapping() const
    { return overlapping_; }

    GridIndexType dofIndex() const
    { return dofIdx_; }

    LocalIndexType indexInElement() const
    { return indexInElement_; }

    GridIndexType elementIndex() const
    { return eIdx_; }

    LocalIndexType localDofIndex() const
    { return indexInElement_; }

private:
    GlobalPosition center_;
    GlobalPosition dofPosition_;
    Scalar volume_;
    LocalIndexType indexInElement_;
    GridIndexType eIdx_;
    GridIndexType dofIdx_;
    bool overlapping_;
};

} // end namespace Dumux

#endif
