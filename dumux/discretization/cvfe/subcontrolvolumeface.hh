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
 * \ingroup CvfeDiscretization
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CVFE_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux {

/*!
 * \ingroup CvfeDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the cvfe scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct CvfeDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = CvfeMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup CvfeDiscretization
 * \brief Class for a sub control volume face in the cvfe method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = CvfeDefaultScvfGeometryTraits<GV> >
class CvfeSubControlVolumeFace
: public SubControlVolumeFaceBase<CvfeSubControlVolumeFace<GV, T>, T>
{
    using ThisType = CvfeSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
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
    CvfeSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class Element>
    CvfeSubControlVolumeFace(CornerStorage corners,
                             GlobalPosition normal,
                             const Element& element,
                             const typename Element::Geometry& elemGeometry,
                             GridIndexType scvfIndex,
                             std::vector<LocalIndexType>&& scvIndices,
                             Dune::GeometryType geomType,
                             bool boundary = false)
    : corners_(corners),
      unitOuterNormal_(normal),
      scvfIndex_(scvfIndex),
      scvIndices_(std::move(scvIndices)),
      boundary_(boundary),
      boundaryFlag_{},
      geometry_(std::make_unique<Geometry>(geomType, corners))
    {
        center_ = geometry_->center();
    }

    //! Constructor for boundary scvfs
    template<class Intersection>
    CvfeSubControlVolumeFace(CornerStorage corners,
                             const Intersection& intersection,
                             const typename Intersection::Geometry& isGeometry,
                             LocalIndexType indexInIntersection,
                             GridIndexType scvfIndex,
                             std::vector<LocalIndexType>&& scvIndices,
                             Dune::GeometryType geomType,
                             bool boundary = false)
    : corners_(corners),
      unitOuterNormal_(intersection.centerUnitOuterNormal()),
      scvfIndex_(scvfIndex),
      scvIndices_(std:: move(scvIndices)),
      boundary_(boundary),
      boundaryFlag_{intersection},
      geometry_(std::make_unique<Geometry>(geomType, corners)), center_(0.0)
    {
        center_ = geometry_->center();
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
        return geometry_->volume();
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

    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    const Geometry& geometry() const
    {
        return *geometry_;
    }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    {
        return boundaryFlag_.get();
    }

private:
    CornerStorage corners_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
    std::unique_ptr<Geometry> geometry_;
    GlobalPosition center_;
};

} // end namespace Dumux

#endif
