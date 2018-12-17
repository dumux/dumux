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
 * \ingroup BoxDiscretization
 * \brief the sub control volume for the box scheme
 */
#ifndef DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_BOX_SUBCONTROLVOLUME_HH

#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
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

    template <class ct>
    struct ScvMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim) corners (1<<(dim))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(dim)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, ScvMLGTraits<Scalar>>;
    using CornerStorage = typename ScvMLGTraits<Scalar>::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
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
    using CornerStorage = typename T::CornerStorage;
    enum { dim = Geometry::mydimension };

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxSubControlVolume() = default;

    // the contructor in the box case
    template<class GeometryHelper>
    BoxSubControlVolume(const GeometryHelper& geometryHelper,
                        LocalIndexType scvIdx,
                        GridIndexType elementIndex,
                        GridIndexType dofIndex)
    : corners_(geometryHelper.getScvCorners(scvIdx)),
      center_(0.0),
      volume_(geometryHelper.scvVolume(corners_)),
      elementIndex_(elementIndex),
      localDofIdx_(scvIdx),
      dofIndex_(dofIndex)
    {
        // compute center point
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
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

    //! The geometry of the sub control volume
    // e.g. for integration
    Geometry geometry() const
    {
        return Geometry(Dune::GeometryTypes::cube(dim), corners_);
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
        // The corner list is defined such that the first entry is the vertex itself
        return corners_[0];
    }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    {
        return elementIndex_;
    }

    //! Return the corner for the given local index
    const GlobalPosition& corner(LocalIndexType localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

private:
    CornerStorage corners_;
    GlobalPosition center_;
    Scalar volume_;
    GridIndexType elementIndex_;
    LocalIndexType localDofIdx_;
    GridIndexType dofIndex_;
};

} // end namespace Dumux

#endif
