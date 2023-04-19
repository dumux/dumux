// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief The spatial parameters for a discrete fracture network embedded in an impermeable matrix.
 */

#ifndef DUMUX_TWOP_FRACTURE_TEST_SPATIALPARAMS_HH
#define DUMUX_TWOP_FRACTURE_TEST_SPATIALPARAMS_HH

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/porousmediumflow/2p/model.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial parameters for for a discrete fracture network embedded in an impermeable matrix.
 */
template<class GridGeometry, class Scalar>
class FractureSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           FractureSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                         FractureSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    //! export permeability type
    using PermeabilityType = Scalar;

    FractureSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {}

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-10; }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Specifies how much the domain is extruded at a given position.
     * \param globalPos The global position where to define the extrusion
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    { return 0.1; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The position at which we evaluate
     */
    auto fluidMatrixInteractionAtPos (const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
