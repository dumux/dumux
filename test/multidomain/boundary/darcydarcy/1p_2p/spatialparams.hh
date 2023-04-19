// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The spatial parameters for the incompressible test.
 */

#ifndef DUMUX_MULTIDOMAIN_1P_2P_TEST_SPATIAL_PARAMS_HH
#define DUMUX_MULTIDOMAIN_1P_2P_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model.
 */
template<class GridGeometry, class Scalar>
class TestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TestSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenNoReg<Scalar>;

public:
    using PermeabilityType = Scalar;

    TestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {}

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return 1e-12;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 0.3;
    }

    /*!
     * \brief Defines the temperature at the given position \f$\mathrm{[K]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 283.15;
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

private:
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
