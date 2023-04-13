// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the dam break problem.
 */
#ifndef DUMUX_DAM_BREAK_SPATIAL_PARAMETERS_HH
#define DUMUX_DAM_BREAK_SPATIAL_PARAMETERS_HH

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the dam break test.
 *
 */
template<class GridGeometry, class Scalar>
class DamBreakSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, DamBreakSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = DamBreakSpatialParams<GridGeometry, Scalar>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    DamBreakSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }


    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \return the bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return 0.0;
    }

private:
    static constexpr Scalar gravity_ = 9.81;
};

} // end namespace Dumux

#endif
