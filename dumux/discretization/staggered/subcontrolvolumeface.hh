// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for a staggered grid geometry helper
 */
template<class GridView>
class BaseStaggeredGeometryHelper
{
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    static constexpr int codimIntersection =  1;

public:

    BaseStaggeredGeometryHelper(const Element& element, const GridView& gridView)
    : element_(element)
    , gridView_(gridView)
    { }

    /*!
    * \brief Updates the current face, i.e. sets the correct intersection
    */
    template<class IntersectionMapper>
    void updateLocalFace(const IntersectionMapper& intersectionMapper, const Intersection& intersection)
    {
        intersection_ = intersection;
    }

    /*!
    * \brief Returns the global dofIdx of the intersection itself
    */
   int dofIndex() const
   {
       //TODO: use proper intersection mapper!
       const auto inIdx = intersection_.indexInInside();
       return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, codimIntersection);
   }

   /*!
   * \brief Returns the local index of the face (i.e. the intersection)
   */
   int localFaceIndex() const
   {
       return intersection_.indexInInside();
   }

private:
   Intersection intersection_; //!< The intersection of interest
   const Element element_; //!< The respective element
   const GridView gridView_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the staggered finite volume scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct StaggeredDefaultScvfGeometryTraits
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::Grid::dimension;
    static constexpr int dimWorld = GridView::Grid::dimensionworld;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim-1, dimWorld>;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for a sub control volume face in the staggered method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on.
 */
template<class GV,
         class T = StaggeredDefaultScvfGeometryTraits<GV> >
class StaggeredSubControlVolumeFace
: public SubControlVolumeFaceBase<StaggeredSubControlVolumeFace<GV, T>, T>
{
    using ThisType = StaggeredSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;

    using Scalar = typename T::Scalar;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

public:
    using GlobalPosition = typename T::GlobalPosition;

    //! state the traits public and thus export all types
    using Traits = T;

    // the default constructor
    StaggeredSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection, class GeometryHelper>
    StaggeredSubControlVolumeFace(const Intersection& is,
                                  const typename Intersection::Geometry& isGeometry,
                                  GridIndexType scvfIndex,
                                  const std::vector<GridIndexType>& scvIndices,
                                  const GeometryHelper& geometryHelper)
    : ParentType()
    , area_(isGeometry.volume())
    , center_(isGeometry.center())
    , unitOuterNormal_(is.centerUnitOuterNormal())
    , scvfIndex_(scvfIndex)
    , scvIndices_(scvIndices)
    , boundary_(is.boundary())
    {
        dofIdx_ = geometryHelper.dofIndex();
        localFaceIdx_ = geometryHelper.localFaceIndex();
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The position of the dof living on the face
    const GlobalPosition& dofPosition() const
    {
        return center_;
    }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    {
        // Return center for now
        return center_;
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return area_;
    }

    //! Returns boolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    //! The unit outer normal vector
    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! Index of the inside sub control volume
    GridIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! Index of the outside sub control volume
    GridIndexType outsideScvIdx() const
    {
        return scvIndices_[1];
    }

    //! The global index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    //! The global index of the dof living on this face
    GridIndexType dofIndex() const
    {
        return dofIdx_;
    }

    //! The local index of this sub control volume face
    LocalIndexType localFaceIdx() const
    {
        return localFaceIdx_;
    }

private:
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    std::vector<GridIndexType> scvIndices_;
    bool boundary_;

    GridIndexType dofIdx_;
    LocalIndexType localFaceIdx_;
};

} // end namespace Dumux

#endif
