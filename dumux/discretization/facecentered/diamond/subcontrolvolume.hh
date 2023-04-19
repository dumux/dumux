// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondSubControlVolume
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUME_HH

#include <array>
#include <utility>
#include <typeinfo>

#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/indextraits.hh>
#include "geometryhelper.hh"

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct FaceCenteredDiamondScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, FCDiamondMLGeometryTraits<Scalar>>;

    using CornerStorage = typename FCDiamondMLGeometryTraits<Scalar>::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;

    static constexpr Dune::GeometryType geometryType(Dune::GeometryType elementType)
    {
        if (elementType == Dune::GeometryTypes::hexahedron)
            return Dune::GeometryTypes::pyramid;
        else
            return Dune::GeometryTypes::simplex(dim);
    }
};

/*!
 * \ingroup DiamondDiscretization
 * \brief Face centered diamond subcontrolvolume face
 */
template<class GridView, class T = FaceCenteredDiamondScvGeometryTraits<GridView>>
class FaceCenteredDiamondSubControlVolume
{
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;

public:
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    FaceCenteredDiamondSubControlVolume() = default;

    FaceCenteredDiamondSubControlVolume(const Scalar& volume,
                                        const GlobalPosition& dofPosition,
                                        const GlobalPosition& center,
                                        const LocalIndexType indexInElement,
                                        const GridIndexType eIdx,
                                        const GridIndexType dofIdx)
    : center_(center)
    , dofPosition_(dofPosition)
    , volume_(volume)
    , indexInElement_(indexInElement)
    , eIdx_(eIdx)
    , dofIdx_(dofIdx)
    {}

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The position of the degree of freedom
    const GlobalPosition& dofPosition() const
    { return dofPosition_; }

    Scalar volume() const
    { return volume_; }

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
};

} // end namespace Dumux

#endif
