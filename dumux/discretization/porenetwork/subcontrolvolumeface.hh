// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_PNM_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_PNM_SUBCONTROLVOLUMEFACE_HH

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Default traits class
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct PNMDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

    using Scalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CornerStorage = std::array<Dune::FieldVector<Scalar, dimWorld>, 1>;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim-1, dimWorld>;

};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Class for a sub control volume face for porenetworks
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = PNMDefaultScvfGeometryTraits<GV> >
class PNMSubControlVolumeFace
: public Dumux::SubControlVolumeFaceBase<PNMSubControlVolumeFace<GV, T>, T>
{
    using ThisType = PNMSubControlVolumeFace<GV, T>;
    using ParentType = Dumux::SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    PNMSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    PNMSubControlVolumeFace(const GlobalPosition& center,
                            const GlobalPosition& unitOuterNormal,
                            const Scalar& area,
                            GridIndexType scvfIndex,
                            std::array<LocalIndexType, 2>&& scvIndices)
    : center_(center),
      unitOuterNormal_(unitOuterNormal),
      area_(area),
      scvfIndex_(scvfIndex),
      scvIndices_(std::move(scvIndices))
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

    //! We assume to always have a pore body and not a pore throat at the boundary
    bool boundary() const
    { return false; }

    //! The unit outer normal of the sub control volume face
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! Index of the inside sub control volume for spatial param evaluation
    LocalIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! Index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    LocalIndexType outsideScvIdx() const
    { return scvIndices_[1]; }

    //! The local index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

private:
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    GridIndexType scvfIndex_;
    std::array<LocalIndexType, 2> scvIndices_;
};

} // end namespace Dumux::PoreNetwork

#endif
