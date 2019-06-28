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
 * \ingroup TwoPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */

#ifndef DUMUX_2P_VOLUME_VARIABLES_HH
#define DUMUX_2P_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class Traits>
class TwoPVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPVolumeVariables<Traits> >;
    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;
    using Indices = typename ModelTraits::Indices;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        phase0Idx = FS::phase0Idx,
        phase1Idx = FS::phase1Idx
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();

public:
    //! Export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of fluid state
    using FluidState = typename Traits::FluidState;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const int nPhaseIdx = 1 - wPhaseIdx_;

        mobility_[wPhaseIdx_] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx_))
            / fluidState_.viscosity(wPhaseIdx_);

        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx_))
            / fluidState_.viscosity(nPhaseIdx);

        // porosity calculation over inert volumefraction
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     * \param solidState A container with the current (physical) state of the solid
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities and enthalpies.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = elemSol[scv.localDofIndex()];

        wPhaseIdx_ = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            if (wPhaseIdx_ == phase1Idx)
            {
                fluidState.setSaturation(phase1Idx, priVars[saturationIdx]);
                fluidState.setSaturation(phase0Idx, 1 - priVars[saturationIdx]);
                pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx_));
                fluidState.setPressure(phase1Idx, priVars[pressureIdx] - pc_);
            }
            else
            {
                const auto Sn = Traits::SaturationReconstruction::reconstructSn(problem.spatialParams(), element,
                                                                                scv, elemSol, priVars[saturationIdx]);
                fluidState.setSaturation(phase1Idx, Sn);
                fluidState.setSaturation(phase0Idx, 1 - Sn);
                pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx_));
                fluidState.setPressure(phase1Idx, priVars[pressureIdx] + pc_);
            }
        }
        else if (formulation == TwoPFormulation::p1s0)
        {
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            if (wPhaseIdx_ == phase1Idx)
            {
                const auto Sn = Traits::SaturationReconstruction::reconstructSn(problem.spatialParams(), element,
                                                                                scv, elemSol, priVars[saturationIdx]);
                fluidState.setSaturation(phase0Idx, Sn);
                fluidState.setSaturation(phase1Idx, 1 - Sn);
                pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx_));
                fluidState.setPressure(phase0Idx, priVars[pressureIdx] + pc_);
            }
            else
            {
                fluidState.setSaturation(phase0Idx, priVars[saturationIdx]);
                fluidState.setSaturation(phase1Idx, 1 - priVars[saturationIdx]);
                pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx_));
                fluidState.setPressure(phase0Idx, priVars[pressureIdx] - pc_);
            }
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numFluidPhases(); ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

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
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
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
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
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
    const int wettingPhaseIdx() const
    {  return wPhaseIdx_; }

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    int wPhaseIdx_;
    Scalar pc_;
    Scalar porosity_;
    PermeabilityType permeability_;
    Scalar mobility_[ModelTraits::numFluidPhases()];
};

} // end namespace Dumux

#endif
