// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief the sub control volume for the cvfe scheme
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_SUBCONTROLVOLUME_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the pq1bubble scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PQ1BubbleDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = PQ1BubbleMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief the sub control volume for the pq1bubble scheme
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GridView, class T = PQ1BubbleDefaultScvGeometryTraits<GridView>>
class PQ1BubbleSubControlVolume
{
    using GlobalPosition = typename T::GlobalPosition;
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    PQ1BubbleSubControlVolume() = default;

    PQ1BubbleSubControlVolume(const Scalar& volume,
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
