// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c chemical nonequilibrium problem.
 */

#ifndef DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH
#define DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams.hh>
#include <dumux/porousmediumflow/fvspatialparamsnonequilibrium.hh>

#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/interfacialarea/interfacialarea.hh>
#include <dumux/material/fluidmatrixinteractions/2p/interfacialarea/pcmax.hh>

namespace Dumux {

/**
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c chemical nonequilibrium problem.
 */
template<class GridGeometry, class Scalar>
class TwoPTwoCChemicalNonequilibriumSpatialParams
: public FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, Scalar,
                                       TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ThisType = TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, Scalar, ThisType>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum {dimWorld=GridView::dimensionworld};

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using WettingNonwettingInterfacialArea = FluidMatrix::InterfacialArea<Scalar,
                                                                          FluidMatrix::InterfacialAreaPcMax,
                                                                          FluidMatrix::WettingNonwettingInterfacialAreaPcSw>;
public:
    //! Export permeability type
    using PermeabilityType = Scalar;

    TwoPTwoCChemicalNonequilibriumSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(1e-11)
    , porosity_(0.4)
    , pcKrSwCurve_("SpatialParams")
    {
        characteristicLength_ = getParam<Scalar>("SpatialParams.MeanPoreSize");
        factorMassTransfer_ = getParam<Scalar>("SpatialParams.MassTransferFactor");
        temperature_ = 273.15 + 25; // -> 25°C

        auto anwParams = WettingNonwettingInterfacialArea::makeBasicParams("SpatialParams.WettingNonwettingArea");

        // determine maximum capillary pressure for wetting-nonwetting surface
        anwParams.setPcMax(pcKrSwCurve_.pc(/*sw = */0.0));

        aNw_ = std::make_unique<WettingNonwettingInterfacialArea>(anwParams);
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the sub-control volume.
     * \return the material parameters object
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_, *aNw_);
    }

    /*!
     * \brief Returns the characteristic length for the mass transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    { return characteristicLength_ ; }

    /*!
     * \brief Return the pre factor the the energy transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    { return factorMassTransfer_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    /*!
     * \brief Returns the temperature in the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos (const GlobalPosition& globalPos) const
    {
        return temperature_;
    }

private:

    const Scalar permeability_;
    const Scalar porosity_;
    static constexpr Scalar eps_ = 1e-6;

    const PcKrSwCurve pcKrSwCurve_;
    std::unique_ptr<const WettingNonwettingInterfacialArea> aNw_;

    Scalar factorMassTransfer_ ;
    Scalar characteristicLength_ ;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
