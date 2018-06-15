// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup CCDiscretization
 * \brief Sub control volumes for cell-centered discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH

#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>
#include <dumux/common/optional.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct CCDefaultScvGeometryTraits
{
    using Geometry = typename GridView::template Codim<0>::Geometry;
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = unsigned int;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
};

/*!
 * \ingroup CCDiscretization
 * \brief Sub control volumes for cell-centered discretization schemes
 * \tparam GV the type of the grid view
 * \tparam T the scv geometry traits
 */
template<class GV,
         class T = CCDefaultScvGeometryTraits<GV> >
class CCSubControlVolume
: public SubControlVolumeBase<CCSubControlVolume<GV, T>, T>
{
    using ThisType = CCSubControlVolume<GV, T>;
    using ParentType = SubControlVolumeBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;

    // In the following, the correct parameter type for the geometry passed to
    // the constructor below is determined. It depends upon whether the
    // `geometry()` method of the `Element` returns a copy or a const reference.
    // In the first case, the correct type is `Geometry&&`, while it is
    // `Geometry` for the second case. Although returning by copy is prescribed
    // by the Dune interface, the grid implementation CpGrid uses a const
    // reference as of Opm 2018.04. Once this is fixed, the parameter type can
    // be hardcoded to `Geometry&&`.
    using Element = typename GV::template Codim<0>::Entity;
    using GeometryRT = typename std::result_of<decltype(&Element::geometry)(Element)>::type;
    using GeometryRTWithoutRef = typename std::remove_reference<GeometryRT>::type;
    static constexpr bool grtIsReference = !std::is_same<GeometryRT, GeometryRTWithoutRef>::value;
    static constexpr bool grtIsConst = std::is_const<GeometryRTWithoutRef>::value;
    using GeometryParamType = std::conditional_t<grtIsReference && grtIsConst, Geometry, Geometry&&>;
public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    CCSubControlVolume() = default;

    // See the explanation above for deriving `GeometryParamType`.
    CCSubControlVolume(GeometryParamType geometry,
                       GridIndexType elementIndex)
    : ParentType()
    , geometry_(std::move(geometry))
    , center_(geometry_.value().center())
    , elementIndex_(elementIndex)
    {}

    //! The copy constrcutor
    CCSubControlVolume(const CCSubControlVolume& other) = default;

    //! The move constrcutor
    CCSubControlVolume(CCSubControlVolume&& other) = default;

    //! The copy assignment operator
    CCSubControlVolume& operator=(const CCSubControlVolume& other)
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(other.geometry_.value());
        center_ = other.center_;
        elementIndex_ = other.elementIndex_;
        return *this;
    }

    //! The move assignment operator
    CCSubControlVolume& operator=(CCSubControlVolume&& other) noexcept
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(std::move(other.geometry_.value()));
        center_ = std::move(other.center_);
        elementIndex_ = std::move(other.elementIndex_);
        return *this;
    }

    //! The center of the sub control volume
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return geometry().volume();
    }

    //! The geometry of the sub control volume
    // e.g. for integration
    const Geometry& geometry() const
    {
        assert((geometry_));
        return geometry_.value();
    }

    //! The index of the dof this scv is embedded in (the global index of this scv)
    GridIndexType dofIndex() const
    {
        return elementIndex();
    }

    //! The element-local index of the dof this scv is embedded in
    LocalIndexType localDofIndex() const
    {
        return 0;
    }

    //! The element-local index of this scv.
    //! In cell-centered schemes there is always only one scv per element.
    LocalIndexType indexInElement() const
    {
        return 0;
    }

    // The position of the dof this scv is embedded in
    const GlobalPosition& dofPosition() const
    {
        return center_;
    }

    //! The global index of the element this scv is embedded in
    GridIndexType elementIndex() const
    {
        return elementIndex_;
    }

    //! Return the corner for the given local index
    GlobalPosition corner(LocalIndexType localIdx) const
    {
        assert(localIdx < geometry().corners() && "provided index exceeds the number of corners");
        return geometry().corner(localIdx);
    }

private:
    // Work around the fact that geometry is not default constructible
    Optional<Geometry> geometry_;
    GlobalPosition center_;
    GridIndexType elementIndex_;
};

} // end namespace Dumux

#endif
