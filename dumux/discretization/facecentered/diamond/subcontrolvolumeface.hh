// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_SUBCONTROLVOLUMEFACE_HH

#include <array>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>


#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct FaceCenteredDiamondScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, FCDiamondMLGeometryTraits<Scalar>>;
    using CornerStorage = typename FCDiamondMLGeometryTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;

    static constexpr Dune::GeometryType interiorGeometryType(Dune::GeometryType)
    { return Dune::GeometryTypes::simplex(dim-1); }
};

/*!
 * \ingroup DiamondDiscretization
 * \brief The SCVF implementation for diamond
 */
template<class GridView, class T = FaceCenteredDiamondScvfGeometryTraits<GridView>>
class FaceCenteredDiamondSubControlVolumeFace
{
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    using GlobalPosition = typename T::GlobalPosition;

    FaceCenteredDiamondSubControlVolumeFace() = default;

    // interior scvf
    FaceCenteredDiamondSubControlVolumeFace(const GlobalPosition& center,
                                            const Scalar area,
                                            const GlobalPosition& normal,
                                            const std::array<LocalIndexType, 2>& scvIndices,
                                            const LocalIndexType localScvfIdx)
    : center_(center)
    , unitOuterNormal_(normal)
    , area_(area)
    , localScvfIdx_(localScvfIdx)
    , scvIndices_(scvIndices)
    , boundary_(false)
    , boundaryFlag_{}
    {}

    // boundary scvf
    FaceCenteredDiamondSubControlVolumeFace(const GlobalPosition& center,
                                            const Scalar area,
                                            const GlobalPosition& normal,
                                            const std::array<LocalIndexType, 2>& scvIndices,
                                            const LocalIndexType localScvfIdx,
                                            const BoundaryFlag& bFlag)
    : center_(center)
    , unitOuterNormal_(normal)
    , area_(area)
    , localScvfIdx_(localScvfIdx)
    , scvIndices_(scvIndices)
    , boundary_(true)
    , boundaryFlag_(bFlag)
    {}

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point of the sub control volume face
    const GlobalPosition& ipGlobal() const
    { return center_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    { return unitOuterNormal_; }

    //! Index of the inside sub control volume
    GridIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! index of the outside sub control volume
    GridIndexType outsideScvIdx() const
    { return scvIndices_[1]; }

    std::size_t numOutsideScvs() const
    { return static_cast<std::size_t>(!boundary()); }

    GridIndexType index() const
    { return localScvfIdx_; }

    bool boundary() const
    { return boundary_; }

    Scalar area() const
    { return area_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    {
        return boundaryFlag_.get();
    }

private:
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    LocalIndexType localScvfIdx_;
    std::array<LocalIndexType, 2> scvIndices_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
};


} // end namespace Dumux

#endif
