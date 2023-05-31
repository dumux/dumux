// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
    , area_(isGeometry.volume())
    , center_(isGeometry.center())
    , unitOuterNormal_(is.centerUnitOuterNormal())
    , scvfIndex_(scvfIndex)
    , scvIndices_(scvIndices)
    , boundary_(isBoundary)
    , boundaryFlag_{is}
    {}

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

    //! returns true if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    //! The unit outer normal of the sub control volume face
    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume
    GridIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! Index of the i-th outside sub control volume or boundary scv index.
    // Results in undefined behaviour if i >= numOutsideScvs()
    GridIndexType outsideScvIdx(int i = 0) const
    {
        return scvIndices_[i+1];
    }

    //! The number of scvs on the outside of this face
    std::size_t numOutsideScvs() const
    {
        return scvIndices_.size()-1;
    }

    //! The global index of this sub control volume face
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
