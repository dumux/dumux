// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup ThreePModel
 * \brief Contains the quantities which are constant within a finite volume in the three-phase model.
 */
#ifndef DUMUX_3P_VOLUME_VARIABLES_HH
#define DUMUX_3P_VOLUME_VARIABLES_HH

#include <dumux/common/valgrind.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Contains the quantities which are constant within a finite volume in the three-phase model.
 */
template <class Traits>
class ThreePVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, ThreePVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, ThreePVolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using Indices = typename Traits::ModelTraits::Indices;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numComponents();

    enum {
        numPhases = Traits::ModelTraits::numPhases(),

        wPhaseIdx = FS::wPhaseIdx,
        gPhaseIdx = FS::gPhaseIdx,
        nPhaseIdx = FS::nPhaseIdx,

        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,
        pressureIdx = Indices::pressureIdx
    };

public:
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Update all quantities for a given control volume
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

        // capillary pressure parameters
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

         // mobilities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            mobility_[phaseIdx] = MaterialLaw::kr(materialParams, phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx),
                                 fluidState_.saturation(gPhaseIdx))
                                 / fluidState_.viscosity(phaseIdx);
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // porosity
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        Valgrind::CheckDefined(permeability_);
    }

    /*!
     * \brief Complete the fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The problem
     * \param element The element
     * \param scv The sub control volume
     * \param fluidState The fluid state
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

        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = elemSol[scv.localDofIndex()];

        const Scalar sw = priVars[swIdx];
        const Scalar sn = priVars[snIdx];
        const Scalar sg = 1.0 - sw - sn;

        Valgrind::CheckDefined(sg);

        fluidState.setSaturation(wPhaseIdx, sw);
        fluidState.setSaturation(gPhaseIdx, sg);
        fluidState.setSaturation(nPhaseIdx, sn);

        /* now the pressures */
        const Scalar pg = priVars[pressureIdx];

        // calculate capillary pressures
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const Scalar pcgw = MaterialLaw::pcgw(materialParams, sw);
        const Scalar pcnw = MaterialLaw::pcnw(materialParams, sw);
        const Scalar pcgn = MaterialLaw::pcgn(materialParams, sw + sn);

        const Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn);
        const Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

        const Scalar pn = pg- pcAlpha * pcgn - (1.0 - pcAlpha)*(pcgw - pcNW1);
        const Scalar pw = pn - pcAlpha * pcnw - (1.0 - pcAlpha)*pcNW1;

        fluidState.setPressure(wPhaseIdx, pw);
        fluidState.setPressure(gPhaseIdx, pg);
        fluidState.setPressure(nPhaseIdx, pn);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // compute and set the viscosity
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);

            // compute and set the density
            const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            const Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
            fluidState.setMolarDensity(phaseIdx, rhoMolar);

            // compute and set the enthalpy
            const Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:
    FluidState fluidState_;
    SolidState solidState_;


private:
    PermeabilityType permeability_;
    Scalar mobility_[numPhases];
};

} // end namespace Dumux

#endif
