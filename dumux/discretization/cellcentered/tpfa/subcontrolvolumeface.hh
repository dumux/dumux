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
 * \ingroup CCTpfaDiscretization
 * \brief The sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <vector>

#include <dune/common/reservedvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct CCTpfaDefaultScvfGeometryTraits
{
    using Grid = typename GridView::Grid;

    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    using Scalar = typename Grid::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using GridIndexStorage = typename std::conditional_t< (dim<dimWorld),
                                                          std::vector<GridIndexType>,
                                                          Dune::ReservedVector<GridIndexType, 2> >;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };
    };

    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar> >;
    using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
    using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The sub control volume face
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = CCTpfaDefaultScvfGeometryTraits<GV> >
class CCTpfaSubControlVolumeFace
: public SubControlVolumeFaceBase<CCTpfaSubControlVolumeFace<GV, T>, T>
{
    using ThisType = CCTpfaSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using GridIndexType = typename T::GridIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using GridIndexStorage = typename T::GridIndexStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    // the default constructor
    CCTpfaSubControlVolumeFace() = default;

    /*!
     * \brief Constructor with intersection
     *
     * \param is The intersection
     * \param isGeometry The geometry of the intersection
     * \param scvfIndex The global index of this scv face
     * \param scvIndices The inside/outside scv indices connected to this face
     * \param isBoundary Bool to specify whether or not the scvf is on an interior or the domain boundary
     */
    template <class Intersection>
    CCTpfaSubControlVolumeFace(const Intersection& is,
                               const typename Intersection::Geometry& isGeometry,
                               GridIndexType scvfIndex,
                               const GridIndexStorage& scvIndices,
                               bool isBoundary)
    : ParentType()
    , geomType_(isGeometry.type())
    , area_(isGeometry.volume())
    , center_(isGeometry.center())
    , unitOuterNormal_(is.centerUnitOuterNormal())
    , scvfIndex_(scvfIndex)
    , scvIndices_(scvIndices)
    , boundary_(isBoundary)
    , boundaryFlag_{is}
    {
        corners_.resize(isGeometry.corners());
        for (int i = 0; i < isGeometry.corners(); ++i)
            corners_[i] = isGeometry.corner(i);
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
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

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    //! The unit outer normal of the sub control volume face
    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    GridIndexType outsideScvIdx(int i = 0) const
    {
        return scvIndices_[i+1];
    }

    //! The number of outside scvs connection via this scv face
    std::size_t numOutsideScvs() const
    {
        return scvIndices_.size()-1;
    }

    //! The global index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    //! return the i-th corner of this sub control volume face
    const GlobalPosition& corner(int i) const
    {
        assert(i < corners_.size() && "provided index exceeds the number of corners");
        return corners_[i];
    }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    {
        return Geometry(geomType_, corners_);
    }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    {
        return boundaryFlag_.get();
    }

private:
    Dune::GeometryType geomType_;
    CornerStorage corners_;
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    GridIndexStorage scvIndices_;
    bool boundary_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux

#endif
