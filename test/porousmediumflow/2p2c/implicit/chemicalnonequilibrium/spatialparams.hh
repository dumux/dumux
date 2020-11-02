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
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c chemical nonequilibrium problem.
 */

#ifndef DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH
#define DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/spatialparams/fvnonequilibrium.hh>

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
: public FVNonEquilibriumSpatialParams<GridGeometry, Scalar,
                                       TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVNonEquilibriumSpatialParams<GridGeometry, Scalar,
                                                     TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>>;

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

private:

    const Scalar permeability_;
    const Scalar porosity_;
    static constexpr Scalar eps_ = 1e-6;

    const PcKrSwCurve pcKrSwCurve_;
    std::unique_ptr<const WettingNonwettingInterfacialArea> aNw_;

    Scalar factorMassTransfer_ ;
    Scalar characteristicLength_ ;
};

} // end namespace Dumux

#endif
