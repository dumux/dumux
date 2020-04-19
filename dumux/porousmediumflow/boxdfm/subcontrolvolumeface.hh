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
 * \ingroup BoxDFMModel
 * \brief The sub control volume face class for the box discrete fracture model.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the box discrete fracture scheme
 *
 * \tparam GV the type of the grid view
 *
 * \note We define new traits for the box-dfm sub-control volume face
 *       as we use a different type of container for storing the scvf corners!
 */
template<class GridView>
struct BoxDfmDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        // However, on fracture scvs the number might be smaller (use ReservedVector)
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
    using LocalIndexType = unsigned int;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar>>;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Class for a sub control volume face in the box discrete fracture method, i.e a
 *        part of the boundary of a sub control volume we compute fluxes on.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = BoxDfmDefaultScvfGeometryTraits<GV> >
class BoxDfmSubControlVolumeFace
: public SubControlVolumeFaceBase<BoxDfmSubControlVolumeFace<GV, T>, T>
{
    using ThisType = BoxDfmSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using GlobalPosition = typename T::GlobalPosition;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

    static_assert(T::dim == 2 || T::dim == 3, "Box-Dfm sub-control volume face only implemented in 2d or 3d");

public:
    //! State the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxDfmSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class GeometryHelper, class Element>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Element& element,
                               const typename Element::Geometry& elemGeometry,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices)
    : corners_(geometryHelper.getScvfCorners(scvfIndex))
    , center_(0.0)
    , unitOuterNormal_(geometryHelper.normal(corners_, scvIndices))
    , area_(geometryHelper.scvfArea(corners_))
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(false)
    , isFractureScvf_(false)
    , boundaryFlag_{}
    , facetIdx_(0)
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for boundary scvfs
    template<class GeometryHelper, class Intersection>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Intersection& intersection,
                               const typename Intersection::Geometry& isGeometry,
                               LocalIndexType indexInIntersection,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices)
    : corners_(geometryHelper.getBoundaryScvfCorners(intersection, isGeometry, indexInIntersection))
    , center_(0.0)
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , area_(geometryHelper.scvfArea(corners_))
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(true)
    , isFractureScvf_(false)
    , boundaryFlag_{intersection}
    , facetIdx_(0)
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for inner fracture scvfs
    template<class GeometryHelper, class Intersection>
    BoxDfmSubControlVolumeFace(const GeometryHelper& geometryHelper,
                               const Intersection& intersection,
                               const typename Intersection::Geometry& isGeometry,
                               LocalIndexType indexInIntersection,
                               GridIndexType scvfIndex,
                               std::vector<LocalIndexType>&& scvIndices,
                               bool boundary)
    : corners_(geometryHelper.getFractureScvfCorners(intersection, isGeometry, indexInIntersection))
    , center_(0.0)
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , boundary_(boundary)
    , isFractureScvf_(true)
    , boundaryFlag_{intersection}
    , facetIdx_(intersection.indexInInside())
    {
        // The area here is given in meters. In order to
        // get the right dimensions, the user has to provide
        // the appropriate aperture in the problem (via an extrusion factor)
        if (T::dim == 3)
            area_ = (corners_[1]-corners_[0]).two_norm();
        else if (T::dim == 2)
            area_ = 1.0;

        // obtain the unit normal vector
        unitOuterNormal_ = geometryHelper.fractureNormal(corners_, intersection, indexInIntersection);

        // compute the scvf center
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return center_; }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    { return boundary_; }

    //! returns the unit normal vector pointing outwards
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! The global index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

    //! Return if this is a fracture scvf
    bool isOnFracture() const
    { return isFractureScvf_; }

    //! The element-local facet index for which a fracture scv was created
    LocalIndexType facetIndexInElement() const
    { assert(isFractureScvf_); return facetIdx_; }

    //! Returns the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

    //! Index of the inside sub control volume for spatial param evaluation
    LocalIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! Index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    LocalIndexType outsideScvIdx() const
    {
        assert(!boundary());
        return scvIndices_[1];
    }

    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    {
        if (isFractureScvf_)
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extract the corners from this scv "
                                                    "and extrude them by the corresponding aperture. ");

        return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_);
    }

private:
    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;
    bool boundary_;
    bool isFractureScvf_;
    BoundaryFlag boundaryFlag_;
    LocalIndexType facetIdx_;
};

} // end namespace Dumux

#endif
