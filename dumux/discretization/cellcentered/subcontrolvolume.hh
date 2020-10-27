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
 * \ingroup CCDiscretization
 * \brief Sub control volumes for cell-centered discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH

#include <memory>

#include <dune/common/fvector.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumebase.hh>

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
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
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
    // `geometry()` method of the `Element` returns a copy or a reference.
    // In the first case, the correct type is `Geometry&&`, while it is
    // `Geometry` for the second case. Although returning by copy is prescribed
    // by the Dune interface, the grid implementation CpGrid uses a const
    // reference as of Opm 2018.04. Once this is fixed, the parameter type can
    // be hardcoded to `Geometry&&` again.
    using Element = typename GV::template Codim<0>::Entity;
    using GeometryRT = decltype(std::declval<Element>().geometry());
    static constexpr bool grtIsReference = std::is_lvalue_reference<GeometryRT>::value;
    using GeometryParamType = std::conditional_t<grtIsReference, Geometry, Geometry&&>;
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
    , geometry_(std::make_unique<Geometry>(std::move(geometry)))
    , center_(geometry_->center())
    , elementIndex_(elementIndex)
    {}

    //! The copy constrcutor
    CCSubControlVolume(const CCSubControlVolume& other)
    { deepCopy_(other); }

    //! The move constrcutor
    CCSubControlVolume(CCSubControlVolume&& other) = default;

    //! The copy assignment operator
    CCSubControlVolume& operator=(const CCSubControlVolume& other)
    {
        deepCopy_(other);
        return *this;
    }

    //! The move assignment operator
    CCSubControlVolume& operator=(CCSubControlVolume&& other) = default;

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
        return *geometry_;
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
    void deepCopy_(const CCSubControlVolume& other)
    {
        if (other.geometry_)
            geometry_ = std::make_unique<Geometry>(other.geometry());
        else
            geometry_.reset();
        center_ = other.center_;
        elementIndex_ = other.elementIndex_;
    }

    // Work around the fact that geometry is not default-constructible and not copy-assignable
    std::unique_ptr<Geometry> geometry_;
    GlobalPosition center_;
    GridIndexType elementIndex_;
};

} // end namespace Dumux

#endif
