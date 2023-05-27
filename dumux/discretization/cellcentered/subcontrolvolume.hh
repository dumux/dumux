// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCDiscretization
 * \brief Sub control volumes for cell-centered discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_SUBCONTROLVOLUME_HH

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
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    CCSubControlVolume() = default;

    template<class Geometry>
    CCSubControlVolume(Geometry&& geometry,
                       GridIndexType elementIndex)
    : ParentType()
    , volume_(geometry.volume())
    , center_(geometry.center())
    , elementIndex_(elementIndex)
    {}

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

private:
    Scalar volume_;
    GlobalPosition center_;
    GridIndexType elementIndex_;
};

} // end namespace Dumux

#endif
