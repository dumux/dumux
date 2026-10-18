// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Contains the quantities which are constant within a finite volume in the two-phase VE model.
 */

#ifndef DUMUX_TWOPVE_VOLUMEVARIABLES_HH
#define DUMUX_TWOPVE_VOLUMEVARIABLES_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/2pve/quantityreconstruction.hh>

namespace Dumux {

//! Traits for isothermal VE flow, without solid-state or energy properties.
template<class PV, class FS, class FST, class PT, class MT>
struct TwoPVEVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FS;
    using FluidState = FST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
};

/*!
 * \ingroup TwoPVEModel
 * \brief Contains the quantities which are constant within a coarse-level finite volume in the two-phase VE model
 *
 * The coarse-level mobilities are the permeability-weighted vertical averages of the reconstructed fine-level mobilities.
 * Permeability [m^2] and dimensionless porosity are column averages; extrusion supplies the column height.
 */
template <class Traits>
class TwoPVEVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;
    using Idx = typename ModelTraits::Indices;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;
    static constexpr int pressureIdx = Idx::pressureIdx;
    static constexpr int saturationIdx = Idx::saturationIdx;
    static constexpr int phase0Idx = FS::phase0Idx;
    static constexpr int phase1Idx = FS::phase1Idx;

    static constexpr auto formulation = ModelTraits::priVarFormulation();
    static_assert(formulation == TwoPFormulation::p0s1, "TwoPVEVolumeVariables only supports the p0s1 formulation");
    static_assert(!ModelTraits::enableEnergyBalance(), "TwoPVEVolumeVariables only supports isothermal flow");

    using GasPlumeDistances = TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = TwoPVE::PhaseViscositiesData<Scalar>;
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
        ParentType::update(elemSol, problem, element, scv);

        const auto columnIdx = problem.gridGeometry().elementMapper().index(element);
        const auto& fineLevelView = problem.fineLevelView();
        const auto columnState = fineLevelView.makeColumnState(element, this->priVars(), problem.spatialParams());
        const Scalar deltaZ = fineLevelView.fineCellHeight();
        constexpr int dim = std::decay_t<decltype(fineLevelView.gridGeometry())>::GridView::dimensionworld;

        const auto& column = fineLevelView.columnMap().column(columnIdx);

        // the wetting-phase pressure determines all coarse-level densities and viscosities
        fluidState_.setDensity(phase0Idx, columnState.densityW);
        fluidState_.setViscosity(phase0Idx, columnState.viscosityW);
        fluidState_.setDensity(phase1Idx, columnState.densityNw);
        fluidState_.setViscosity(phase1Idx, columnState.viscosityNw);

        const Scalar pcCoarse = fineLevelView.quantityReconstructor().computeCapillaryPressureCoarse(
             columnState.gasPlumeDistance,
             PhaseDensities{columnState.densityW, columnState.densityNw},
             columnState.gravityNorm,
             columnState.entryPressure);
        const Scalar permeabilityIntegral = problem.spatialParams().permeabilityIntegralAtElement(element);
        columnHeight_ = problem.spatialParams().columnHeightAtElement(element);
        porosity_ = problem.spatialParams().porosityIntegralAtElement(element)/columnHeight_;

        // column integrals of the fine-level permeability times the reconstructed fine-level mobilities
        mobility_[phase0Idx] = 0.0;
        mobility_[phase1Idx] = 0.0;
        const auto& spatialParamsFine = problem.spatialParams().spatialParamsFine();
        for (const auto& fineElement : column)
        {
            // relative height instead of absolute height is required for comparison with gas plume distance
            const Scalar fineElementHeight = fineElement.geometry().center()[dim - 1] - fineLevelView.gridGeometry().bBoxMin()[dim - 1];

            const auto reconstructedMobilities = fineLevelView.quantityReconstructor().reconstructMobilitiesFine(
                 GasPlumeDistances{columnState.gasPlumeDistance, columnState.minimumGasPlumeDistance},
                 PhaseDensities{columnState.densityW, columnState.densityNw},
                 PhaseViscosities{columnState.viscosityW, columnState.viscosityNw},
                 columnState.gravityNorm,
                 fineElementHeight,
                 deltaZ,
                 BrooksCoreyParameters{columnState.brooksCoreyLambda,columnState.entryPressure});

            const Scalar permeabilityFine = spatialParamsFine.permeabilityAtElement(fineElement);
            mobility_[phase0Idx] += permeabilityFine * reconstructedMobilities[phase0Idx] * deltaZ;
            mobility_[phase1Idx] += permeabilityFine * reconstructedMobilities[phase1Idx] * deltaZ;
        }

        completeFluidStateCoarse(elemSol, problem, element, scv, fluidState_, columnState.temperature, pcCoarse, columnState.gasPlumeDistance);

        // permeability-weighted column average, the denominator is the column integral of the fine-level permeability
        mobility_[phase0Idx] /= permeabilityIntegral;
        mobility_[phase1Idx] /= permeabilityIntegral;

        permeability_ = permeabilityIntegral/columnHeight_;
    }

    /*!
     * \brief Sets complete fluid state. This function is used for the coarse-level elements.
     *
     * \param elemSol       a vector containing all primary variables connected to the element
     * \param problemCoarse the object specifying the coarse-level problem which ought to be simulated
     * \param element       an element which contains part of the control volume
     * \param scv           the sub-control volume
     * \param fluidState    a container with the current (physical) state of the fluid
     * \param temperature   prescribed column temperature
     * \param pcCoarse      coarsened capillary pressure of coarse-level element
     * \param gasPlumeDist  gas plume distance/height belonging to the coarse-level element
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidStateCoarse(const ElemSol& elemSol,
                                  const Problem& problemCoarse,
                                  const Element& element,
                                  const Scv& scv,
                                  FluidState& fluidState,
                                  Scalar temperature,
                                  Scalar pcCoarse,
                                  Scalar gasPlumeDist)
    {
        for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            fluidState.setTemperature(phaseIdx, temperature);
        const auto& priVars = elemSol[scv.localDofIndex()];

        const auto wPhaseIdx = problemCoarse.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        if (wPhaseIdx != phase0Idx)
            DUNE_THROW(Dune::InvalidStateException, "The two-phase VE model requires the first phase to be the wetting phase");
        fluidState.setWettingPhase(wPhaseIdx);

        fluidState.setSaturation(phase1Idx, priVars[saturationIdx]);
        fluidState.setSaturation(phase0Idx, 1.0 - priVars[saturationIdx]);

        pc_ = pcCoarse;
        fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
        fluidState.setPressure(phase1Idx, priVars[pressureIdx] + pc_);

        gasPlumeDist_ = gasPlumeDist;
    }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

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
     * \brief Returns the dimensionless column-average porosity.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the column-average intrinsic permeability in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    //! Vertical integral of porosity [m].
    Scalar porosityIntegral() const
    { return porosity_*columnHeight_; }

    //! Vertical integral of permeability [m^3].
    PermeabilityType permeabilityIntegral() const
    { return permeability_*columnHeight_; }

    /*!
     * \brief Returns the wetting phase index
     */
    int wettingPhase() const
    { return fluidState_.wettingPhase(); }

    /*!
     * \brief Returns the gas plume distance within a coarse-level element
     */
    Scalar gasPlumeDist() const
    { return gasPlumeDist_; }

    //! Height of the reconstructed column [m].
    Scalar columnHeight() const { return columnHeight_; }

protected:
    FluidState fluidState_;

private:
    Scalar pc_;
    PermeabilityType permeability_;
    Scalar porosity_;
    Scalar columnHeight_;
    Scalar mobility_[ModelTraits::numFluidPhases()];

    Scalar gasPlumeDist_;
};

} // end namespace Dumux

#endif
