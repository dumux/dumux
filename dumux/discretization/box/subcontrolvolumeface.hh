// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/geometry/center.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the box scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct BoxDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = BoxMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDefaultScvfGeometryTraits<GV> >
class BoxSubControlVolumeFace
: public SubControlVolumeFaceBase<BoxSubControlVolumeFace<GV, T>, T>
{
    using ThisType = BoxSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;
    static constexpr int dim = Geometry::mydimension;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class Corners, class Element>
    BoxSubControlVolumeFace(const Corners& corners,
                            const GlobalPosition& normal,
                            const Element& element,
                            const typename Element::Geometry& elemGeometry,
                            GridIndexType scvfIndex,
                            std::vector<LocalIndexType>&& scvIndices,
                            bool boundary = false)
    : center_(Dumux::center(corners))
    , unitOuterNormal_(normal)
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(boundary)
    , boundaryFlag_{}
    {
        area_ = Dumux::convexPolytopeVolume<dim>(
            Dune::GeometryTypes::cube(dim),
            [&](unsigned int i){ return corners[i]; }
        );
    }

    //! Constructor for boundary scvfs
    template<class Corners, class Intersection>
    BoxSubControlVolumeFace(const Corners& corners,
                            const GlobalPosition& normal,
                            const Intersection& intersection,
                            const typename Intersection::Geometry& isGeometry,
                            LocalIndexType indexInIntersection,
                            GridIndexType scvfIndex,
                            std::vector<LocalIndexType>&& scvIndices,
                            bool boundary = false)
    : center_(Dumux::center(corners))
    , unitOuterNormal_(normal)
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(boundary)
    , boundaryFlag_{intersection}
    {
        area_ = Dumux::convexPolytopeVolume<dim>(
            Dune::GeometryTypes::cube(dim),
            [&](unsigned int i){ return corners[i]; }
        );
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    {
        return center_;
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return area_;
    }

    //! returns true if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume
    LocalIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

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

    //! The local index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    {
        return boundaryFlag_.get();
    }

private:
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux

#endif
