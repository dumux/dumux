// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief the sub control volume for the box scheme
 */
#ifndef DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/geometry/center.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the box scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct BoxDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = BoxMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
};

/*!
 * \ingroup BoxDiscretization
 * \brief the sub control volume for the box scheme
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDefaultScvGeometryTraits<GV> >
class BoxSubControlVolume
: public SubControlVolumeBase<BoxSubControlVolume<GV, T>, T>
{
    using ThisType = BoxSubControlVolume<GV, T>;
    using ParentType = SubControlVolumeBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    static constexpr int dim = Geometry::mydimension;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxSubControlVolume() = default;

    // the constructor in the box case
    template<class Corners>
    BoxSubControlVolume(const Corners& corners,
                        LocalIndexType scvIdx,
                        GridIndexType elementIndex,
                        GridIndexType dofIndex)
    : dofPosition_(corners[0])
    , center_(Dumux::center(corners))
    , elementIndex_(elementIndex)
    , localDofIdx_(scvIdx)
    , dofIndex_(dofIndex)
    {
        // The corner list is defined such that the first entry is the vertex itself
        volume_ = Dumux::convexPolytopeVolume<dim>(
            Dune::GeometryTypes::cube(dim),
            [&](unsigned int i){ return corners[i]; }
        );
    }

    //! The center of the sub control volume
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return volume_;
    }

    //! The element-local index of the dof this scv is embedded in
    LocalIndexType localDofIndex() const
    {
        return localDofIdx_;
    }

    //! The element-local index of this scv.
    //! For the standard box scheme this is the local dof index.
    LocalIndexType indexInElement() const
    {
        return localDofIdx_;
    }

    //! The index of the dof this scv is embedded in
    GridIndexType dofIndex() const
    {
        return dofIndex_;
    }

    // The position of the dof this scv is embedded in
    const GlobalPosition& dofPosition() const
    {
        return dofPosition_;
    }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    {
        return elementIndex_;
    }

private:
    GlobalPosition dofPosition_;
    GlobalPosition center_;
    Scalar volume_;
    GridIndexType elementIndex_;
    LocalIndexType localDofIdx_;
    GridIndexType dofIndex_;
};

} // end namespace Dumux

#endif
