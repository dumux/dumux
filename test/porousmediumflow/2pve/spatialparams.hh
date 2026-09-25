// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief The spatial params for the vertical equilibrium Darcy test.
 */

#ifndef DUMUX_TEST_TWOPVE_SPATIAL_PARAMS_HH
#define DUMUX_TEST_TWOPVE_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>
#include <dumux/porousmediumflow/2pve/spatialparams.hh>

#include "spatialparams_fine.hh"

namespace Dumux {

/*!
 * \ingroup TwoPVETests
 * \brief The coarse-level spatial params for the vertical equilibrium Darcy test.
 */
template<class GridGeometry, class Scalar>
class TwoPTestSpatialParams
: public TwoPVESpatialParams<GridGeometry, Scalar, TwoPTestFineSpatialParams<GridGeometry, Scalar>, TwoPTestSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPTestSpatialParams<GridGeometry, Scalar>;
    using SpatialParamsFine = TwoPTestFineSpatialParams<GridGeometry, Scalar>;
    using ParentType = TwoPVESpatialParams<GridGeometry, Scalar, SpatialParamsFine, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    TwoPTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                          const VEColumnMapping<GridGeometry, Scalar>& columnMapping,
                          std::shared_ptr<const SpatialParamsFine> spatialParamsFine,
                          const Scalar fineCellHeight)
    : ParentType(gridGeometry, columnMapping, spatialParamsFine, fineCellHeight)
    , pcKrSwCurve_("SpatialParams")
    {}

    /*!
     * \brief Returns the parameter object for the capillary-pressure/saturation material law
     *
     * \param globalPos the coordinates
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * \param globalPos the coordinates
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 326.0; // 53°C
    }

private:
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
