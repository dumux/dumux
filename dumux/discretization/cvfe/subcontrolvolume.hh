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
 * \brief the sub control volume for the cvfe scheme
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_CVFE_SUBCONTROLVOLUME_HH

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>

namespace Dumux {

/*!
 * \ingroup CvfeDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cvfe scheme
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct CvfeDefaultScvGeometryTraits
{
    using Grid = typename GridView::Grid;

    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename Grid::ctype;
    using GeometryTraits = CvfeMLGeometryTraits<Scalar>;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, GeometryTraits>;
    using CornerStorage = typename GeometryTraits::template CornerStorage<dim, dimWorld>::Type;
    using GlobalPosition = typename CornerStorage::value_type;
};

/*!
 * \ingroup CvfeDiscretization
 * \brief the sub control volume for the cvfe scheme
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV,
         class T = CvfeDefaultScvGeometryTraits<GV> >
class CvfeSubControlVolume
: public SubControlVolumeBase<CvfeSubControlVolume<GV, T>, T>
{
    using ThisType = CvfeSubControlVolume<GV, T>;
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
    CvfeSubControlVolume() = default;

    // the contructor in the cvfe case
    CvfeSubControlVolume(CornerStorage corners,
                         LocalIndexType scvIdx,
                         GridIndexType elementIndex,
                         GridIndexType dofIndex,
                         GlobalPosition dofPosition,
                         Dune::GeometryType geomType)
    : corners_(corners),
      elementIndex_(elementIndex),
      localDofIdx_(scvIdx),
      dofIndex_(dofIndex),
      dofPosition_(dofPosition),
      geometry_(std::make_unique<Geometry>(geomType, corners))
    { }

    //! The copy constrcutor
    CvfeSubControlVolume(const CvfeSubControlVolume& other)
    { deepCopy_(other); }

    //! The move constrcutor
    CvfeSubControlVolume(CvfeSubControlVolume&& other) = default;

    //! The copy assignment operator
    CvfeSubControlVolume& operator=(const CvfeSubControlVolume& other)
    {
        deepCopy_(other);
        return *this;
    }

    //! The move assignment operator
    CvfeSubControlVolume& operator=(CvfeSubControlVolume&& other) = default;

    //! The center of the sub control volume
    GlobalPosition center() const
    {
        return geometry_->center();
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return geometry_->volume();
    }

    //! The geometry of the sub control volume
    // e.g. for integration
    const Geometry& geometry() const
    {
        return *geometry_;
    }

    //! The element-local index of the dof this scv is embedded in
    LocalIndexType localDofIndex() const
    {
        return localDofIdx_;
    }

    //! The element-local index of this scv.
    //! For the standard cvfe scheme this is the local dof index.
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
        return dofPosition_;
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
    void deepCopy_(const CvfeSubControlVolume& other)
    {
        corners_ = other.corners_;
        localDofIdx_ = other.localDofIdx_;
        elementIndex_ = other.elementIndex_;
        dofIndex_ = other.dofIndex_;
        dofPosition_ = other.dofPosition_;
        if (other.geometry_)
            geometry_ = std::make_unique<Geometry>(other.geometry_->type(), corners_);
        else
            geometry_.reset();
    }

    CornerStorage corners_;
    GridIndexType elementIndex_;
    LocalIndexType localDofIdx_;
    GridIndexType dofIndex_;
    GlobalPosition dofPosition_;
    std::unique_ptr<Geometry> geometry_;
};

} // end namespace Dumux

#endif
