// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Implementation of a boundary face related to primary grid elements (dune intersections)
 */
#ifndef DUMUX_DISCRETIZATION_BOUNDARY_FACE_HH
#define DUMUX_DISCRETIZATION_BOUNDARY_FACE_HH

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>

namespace Dumux::Experimental {

//! Traits for an efficient corner storage
template <class ct>
struct BoundaryFaceMLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the maximum number of corners in advance
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<mydim)>;
    };
};

/*!
 * \ingroup Discretization
 * \brief Default traits class to be used for the boundary faces
 *
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct DefaultBoundaryFaceGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = BoundaryFaceMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup Discretization
 * \brief Class for a boundary face related to primary grid elements (dune intersections)
 *
 * \tparam GV the type of the grid view
 * \tparam T the face geometry traits
 */
template<class GV,
         class T = DefaultBoundaryFaceGeometryTraits<GV> >
class BoundaryFace
{
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
    BoundaryFace() = default;

    //! Constructor for boundary faces
    BoundaryFace(const GlobalPosition& center,
                 const Scalar area,
                 const GlobalPosition& normal,
                 const LocalIndexType localIdx,
                 const LocalIndexType intersectionIdx,
                 const BoundaryFlag& bFlag)
    : center_(center)
    , unitOuterNormal_(normal)
    , area_(area)
    , localIdx_(localIdx)
    , intersectionIdx_(intersectionIdx)
    , boundary_(true)
    , boundaryFlag_(bFlag)
    {}

    //! The center of the face
    const GlobalPosition& center() const
    { return center_; }

    //! The area of the face
    Scalar area() const
    { return area_; }

    bool boundary() const
    { return boundary_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    { return unitOuterNormal_; }

    //! The local index of this face
    LocalIndexType index() const
    { return localIdx_; }

    //! The index of the intersection this face corresponds to (intersection.indexInInside())
    LocalIndexType intersectionIndex() const
    { return intersectionIdx_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

private:
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    LocalIndexType localIdx_;
    LocalIndexType intersectionIdx_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux::Experimental

#endif
