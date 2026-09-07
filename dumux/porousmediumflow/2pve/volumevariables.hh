// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief Contains the quantities which are constant within a finite volume in the two-phase VE model.
 */

#ifndef DUMUX_TWOPVE_VOLUMEVARIABLES_HH
#define DUMUX_TWOPVE_VOLUMEVARIABLES_HH

#include <cstddef>
#include <vector>

#include <dumux/common/properties.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/2pve/quantityreconstruction.hh>

namespace Dumux {

template <class Traits>
class TwoPVEVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPVEVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPVEVolumeVariables<Traits> >;
    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;
    using Idx = typename ModelTraits::Indices;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    enum
    {
        pressureIdx = Idx::pressureIdx,
        saturationIdx = Idx::saturationIdx,

        phase0Idx = FS::phase0Idx,
        phase1Idx = FS::phase1Idx,
        numPhases = FS::numPhases
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();
    static_assert(formulation == TwoPFormulation::p0s1, "TwoPVEVolumeVariables only supports the p0s1 formulation");

    using GasPlumeDistances = TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = TwoPVE::PhaseViscositiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:
    //! Export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! Export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of fluid state
    using FluidState = typename Traits::FluidState;
    //! Export the indices
    using Indices = typename ModelTraits::Indices;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Updates all quantities for a given control volume, these volume variables refer to the coarse level of the VE scheme.
     *
     * \param elemSol a vector containing all primary variables connected to the element
     * \param problem the object specifying the problem which ought to be simulated
     * \param element an element which contains part of the control volume
     * \param scv     the sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
        priVars_ = elemSol[scv.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(element, scv, elemSol);
        const int columnIdx = problem.gridGeometry().elementMapper().index(element);
        // compute column state
        const auto columnState = problem.getFineLevelView()->makeColumnState(element, priVars_, problem.spatialParams());
        const Scalar deltaZ = problem.getFineLevelView()->fineCellHeight();
        unsigned int dim = GlobalPosition::dimension;

        const auto& column = problem.getFineLevelView()->columnMap().column(columnIdx);
        std::vector<Scalar> mobilitiesCoarse(numPhases, 0.0);
        std::vector<Scalar> mobWFineEntries(column.size(), 0.0);
        std::vector<Scalar> mobNwFineEntries(column.size(), 0.0);

        //use wetting-phase pressure for computation of all coarse-level densities and viscosities
        fluidState_.setDensity(phase0Idx, columnState.densityW);
        fluidState_.setViscosity(phase0Idx, columnState.viscosityW);
        fluidState_.setDensity(phase1Idx, columnState.densityNw);
        fluidState_.setViscosity(phase1Idx, columnState.viscosityNw);

        //compute coarse-level capillary pressure
        const Scalar pcCoarse = problem.getFineLevelView()->quantityReconstructor().computeCapillaryPressureCoarse(
             columnState.gasPlumeDistance,
             PhaseDensities{columnState.densityW, columnState.densityNw},
             columnState.gravityNorm,
             columnState.entryPressure);
        const Scalar permeabilityCoarse = problem.spatialParams().permeability(element, scv, elemSol);

        // const auto firstCellIterator = column.cbegin();
        Dumux::parallelFor(column.size(), [&](const std::size_t columnElementIdx)
        {
            const auto& fineElement = column[columnElementIdx];
            const Scalar fineElementHeight = fineElement.geometry().center()[dim - 1] - problem.getFineLevelView()->gridGeometry().bBoxMin()[dim - 1]; // relative height instead of absolute height is required for comparison with gas plume distance

            //calculate fine-level mobilities
            std::vector<Scalar> reconstructedMobilites = problem.getFineLevelView()->quantityReconstructor().reconstMobilitiesFine(
                 GasPlumeDistances{columnState.gasPlumeDistance, columnState.minimumGasPlumeDistance},
                 PhaseDensities{columnState.densityW, columnState.densityNw},
                 PhaseViscosities{columnState.viscosityW, columnState.viscosityNw},
                 ResidualSaturations{columnState.swr,columnState.snr},
                 columnState.gravityNorm,
                 fineElementHeight,
                 deltaZ,
                 BrooksCoreyParameters{columnState.brooksCoreyLambda,columnState.entryPressure});

            mobWFineEntries[columnElementIdx] = problem.spatialParams().spatialParamsFine().permeabilityAtElement(fineElement) * reconstructedMobilites[phase0Idx] * deltaZ;
            mobNwFineEntries[columnElementIdx] = problem.spatialParams().spatialParamsFine().permeabilityAtElement(fineElement) * reconstructedMobilites[phase1Idx] * deltaZ;
        });

        for(int columnElementsIdx=0; columnElementsIdx<column.size(); columnElementsIdx++)
        {
            mobilitiesCoarse[phase0Idx] += mobWFineEntries[columnElementsIdx];
            mobilitiesCoarse[phase1Idx] += mobNwFineEntries[columnElementsIdx];
        }

        completeFluidStateCoarse(elemSol, problem, element, scv, fluidState_, solidState_, pcCoarse, columnState.gasPlumeDistance);

        mobilitiesCoarse[phase0Idx] /= permeabilityCoarse;
        mobilitiesCoarse[phase1Idx] /= permeabilityCoarse;

        mobility_[phase0Idx] = mobilitiesCoarse[phase0Idx];
        mobility_[phase1Idx] = mobilitiesCoarse[phase1Idx];

        // porosity calculation over inert volumefraction
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = permeabilityCoarse;
        EnergyVolVars::updateEffectiveThermalConductivity();
    }


    /*!
     * \brief Sets complete fluid state. This function is used for the coarse-level elements.
     *
     * \param elemSol       a vector containing all primary variables connected to the element
     * \param problemCoarse the object specifying the coarse-level problem which ought to be simulated
     * \param element       an element which contains part of the control volume
     * \param scv           the sub-control volume
     * \param fluidState    a container with the current (physical) state of the fluid
     * \param solidState    a container with the current (physical) state of the solid
     * \param pcCoarse      coarsened capillary pressure of coarse-level element
     * \param gasPlumeDist  gas plume distance/height belonging to the coarse-level element
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities and enthalpies.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidStateCoarse(const ElemSol& elemSol,
                                  const Problem& problemCoarse,
                                  const Element& element,
                                  const Scv& scv,
                                  FluidState& fluidState,
                                  SolidState& solidState,
                                  const Scalar& pcCoarse,
                                  const Scalar& gasPlumeDist)
    {
        EnergyVolVars::updateTemperature(elemSol, problemCoarse, element, scv, fluidState, solidState);
        const auto& spatialParams = problemCoarse.spatialParams();
        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto wPhaseIdx = spatialParams.template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);
        fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
        if (fluidState.wettingPhase() == phase1Idx)
        {
            fluidState.setSaturation(phase1Idx, priVars[saturationIdx]);
            fluidState.setSaturation(phase0Idx, 1 - priVars[saturationIdx]);

            pc_ = pcCoarse;
            fluidState.setPressure(phase1Idx, priVars[pressureIdx] + pc_);

            gasPlumeDist_ = gasPlumeDist;
        }
        else
        {
            const auto Sn = Traits::SaturationReconstruction::reconstructSn(spatialParams, element, scv, elemSol, priVars[saturationIdx]);

            fluidState.setSaturation(phase1Idx, Sn);
            fluidState.setSaturation(phase0Idx, 1 - Sn);

            pc_ = pcCoarse;
            fluidState.setPressure(phase1Idx, priVars[pressureIdx] + pc_);
        }

        gasPlumeDist_ = gasPlumeDist;

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx) {
            // compute and set the enthalpy
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }


    /*!
     * \brief Returns the phase state for the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume
     * in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return pc_; }

    /*!
     * \brief Returns temperature inside the sub-control volume
     * in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx the phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the wetting phase index
     */
    int wettingPhase() const
    {  return fluidState_.wettingPhase(); }

    /*!
     * \brief Returns the gas plume distance within a coarse-level element
     */
    const Scalar& gasPlumedist() const
    { return gasPlumeDist_; }

    /*!
     * \brief Returns how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    Scalar pc_;
    PermeabilityType permeability_;
    Scalar mobility_[ModelTraits::numFluidPhases()];

    Scalar gasPlumeDist_;
    Scalar extrusionFactor_; //extrusionFactor from ParentType is "overloaded"

    PrimaryVariables priVars_;
};

} // end namespace Dumux

#endif
