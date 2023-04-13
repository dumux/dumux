// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowTests
 * \brief Definition of the spatial parameters for the 3D Channel Problems.
 */

#ifndef DUMUX_CHANNEL3D_SPATIAL_PARAMS_HH
#define DUMUX_CHANNEL3D_SPATIAL_PARAMS_HH

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar>
class Channel3DSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, Channel3DSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, Channel3DSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename SubControlVolume::Traits::GlobalPosition;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr bool enablePseudoThreeDWallFriction = dim != 3;

public:
    Channel3DSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {

        height_ = getParam<Scalar>("Problem.Height");

        if constexpr (enablePseudoThreeDWallFriction)
            extrusionFactor_ = 2.0/3.0 * height_;
        else
            extrusionFactor_ = 1.0;

    }

    /*!
     * \brief Return the temperature in the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    { return extrusionFactor_; }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusionFactor_; }

private:
    Scalar extrusionFactor_;
    Scalar height_;
};

} // end namespace Dumux

#endif
