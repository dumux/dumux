// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ???
 * \brief The spatial params the trace operator test
 */

#ifndef DUMUX_TRACE_OPERATOR_SPATIAL_PARAMS_HH
#define DUMUX_TRACE_OPERATOR_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup ???
 * \brief The spatial params the trace operator test
 */
template<class GridGeometry, class Scalar>
class TraceTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             TraceTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using ThisType = TraceTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    TraceTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     * \param globalPos The position inside the domain.
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return 1.0; }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     * \param globalPos The position inside the domain.
     */
    Scalar porosityAtPos(const GlobalPosition&) const
    { return 1.0; }

    /*!
     * \brief Define the temperature in the domain in Kelvin.
     * \param globalPos The position inside the domain.
     */
    Scalar temperatureAtPos(const GlobalPosition&) const
    { return 283.15; }

    /*!
     * \brief Define the extrusion of the domain at a given position.
     * \param globalPos The position inside the domain.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition&) const
    { return 1.0; }
};

} // end namespace Dumux

#endif
