// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief Mimics the volume variables structure of DuMux applied to the fine level of the VE scheme.
 */

#ifndef DUMUX_TWOPVE_FINE_LEVEL_ELEMENTSTATE_HH
#define DUMUX_TWOPVE_FINE_LEVEL_ELEMENTSTATE_HH

#include <array>
#include <vector>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2pve/quantityreconstruction.hh>

namespace Dumux {

template<class Scalar>
struct TwoPVEColumnState
{
    Scalar pwCoarse{};
    Scalar swCoarse{};
    Scalar temperature{};

    Scalar domainHeight{};

    Scalar swr{};
    Scalar snr{};

    Scalar brooksCoreyLambda{};
    Scalar entryPressure{};
    Scalar gravityNorm{};

    Scalar densityW{};
    Scalar densityNw{};
    Scalar viscosityW{};
    Scalar viscosityNw{};

    Scalar gasPlumeDistance{};
    Scalar minimumGasPlumeDistance{}; // required for hysteresis
};



template<class TypeTag>
class TwoPVEFineLevelElementState
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    enum {
        dim = GridView::dimension,
    };
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {
        wettingPhaseIdx = FluidSystem::phase0Idx,
        nonwettingPhaseIdx = FluidSystem::phase1Idx,
        numPhases = FluidSystem::numPhases
    };
    using QuantityReconstructor = TwoPVEQuantityReconst<TypeTag>;

    using GasPlumeDistances = TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = TwoPVE::PhaseViscositiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:

    /*!
     * \brief Updates the state of a fine element
     *
     * \param fineElement       fine-level element
     * \param column            column that contains the fine element
     * \param fineSpatialParams fine-level spatial parameters
     * \param reconstructor     object for reconstructing fine- and coarse-level properties
     * \param fineCellHeight    vertical height (relative to the bottom of the domain) of fine-level element
     */
    template<typename SpatialParamsFine>
    void update(const Element& fineElement,
                const TwoPVEColumnState<Scalar>& column,
                const SpatialParamsFine& fineSpatialParams,
                const QuantityReconstructor& reconstructor,
                const Scalar fineCellHeight)
    {
        const Scalar fineElementPosZ = fineElement.geometry().center()[dim-1];
        const Scalar heightAboveBottom = fineElementPosZ - fineSpatialParams.gridGeometry().bBoxMin()[dim-1]; // relative height instead of absolute height is required for reconstruction functions

        const std::vector<Scalar> pressures = reconstructor.reconstPressure(
            column.gasPlumeDistance,
            PhaseDensities{column.densityW, column.densityNw},
            column.gravityNorm,
            heightAboveBottom,
            column.pwCoarse,
            column.entryPressure);
        pressure_[wettingPhaseIdx] = pressures[wettingPhaseIdx];
        pressure_[nonwettingPhaseIdx] = pressures[nonwettingPhaseIdx];

        const Scalar saturationW = reconstructor.reconstructSaturation(
            GasPlumeDistances{column.gasPlumeDistance, column.minimumGasPlumeDistance},
            PhaseDensities{column.densityW, column.densityNw},
            ResidualSaturations{column.swr, column.snr},
            column.gravityNorm,
            heightAboveBottom,
            fineCellHeight,
            BrooksCoreyParameters{column.brooksCoreyLambda, column.entryPressure});
        saturation_[wettingPhaseIdx] = saturationW;
        saturation_[nonwettingPhaseIdx] = 1.0 - saturationW;

        const Scalar capillaryPressure = reconstructor.reconstCapillaryPressure(
            column.gasPlumeDistance,
            PhaseDensities{column.densityW, column.densityNw},
            column.gravityNorm,
            heightAboveBottom,
            column.entryPressure);
        capillaryPressure_ = capillaryPressure;

        const std::vector<Scalar> mobilites = reconstructor.reconstMobilitiesFine(
            GasPlumeDistances{column.gasPlumeDistance, column.minimumGasPlumeDistance},
            PhaseDensities{column.densityW, column.densityNw},
            PhaseViscosities{column.viscosityW, column.viscosityNw},
            ResidualSaturations{column.swr,column.snr},
            column.gravityNorm,
            heightAboveBottom,
            fineCellHeight,
            BrooksCoreyParameters{column.brooksCoreyLambda, column.entryPressure});
        mobility_[wettingPhaseIdx] = mobilites[wettingPhaseIdx];
        mobility_[nonwettingPhaseIdx] = mobilites[nonwettingPhaseIdx];

        density_[wettingPhaseIdx] = column.densityW;
        density_[nonwettingPhaseIdx] = column.densityNw;

        viscosity_[wettingPhaseIdx] = column.viscosityW;
        viscosity_[nonwettingPhaseIdx] = column.viscosityNw;

        permeability_ = fineSpatialParams.permeabilityAtElement(fineElement);
        porosity_ = fineSpatialParams.porosityAtElement(fineElement);

        temperature_ = column.temperature;
        gasPlumeDistance_ = column.gasPlumeDistance;
    }


    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure within the control volume
     * in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return capillaryPressure_; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    Scalar permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns temperature inside the sub-control volume
     * in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the gas plume distance within a coarse-level element
     */
    Scalar gasPlumeDist() const
    { return gasPlumeDistance_; }

private:
    std::array<Scalar, numPhases> pressure_{};
    std::array<Scalar, numPhases> saturation_{};
    Scalar capillaryPressure_{};
    std::array<Scalar, numPhases> mobility_{};
    std::array<Scalar, numPhases> density_{};
    std::array<Scalar, numPhases> viscosity_{};
    Scalar permeability_{};
    Scalar porosity_{};
    Scalar temperature_{};
    Scalar gasPlumeDistance_{};
};

} // end namespace Dumux

#endif
