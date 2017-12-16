// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 *
 * \brief Contains the volume variables of the kinetic energy transfer
 *        module of the M-phase, N-component model.
 */
#ifndef DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_KINETIC_HH
#define DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_KINETIC_HH

#include <dumux/porousmediumflow/mpnc/implicit/volumevariables.hh>
#include "vtkwriterkinetic.hh"

namespace Dumux
{
/*!
 * \brief Contains the volume variables of the kinetic energy transfer
 *        module of the M-phase, N-component model.
 *        Specialization for the case of *3* energy balance equations.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/3>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { temperature0Idx = Indices::temperature0Idx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { sPhaseIdx = FluidSystem::sPhaseIdx };
    enum { numEnergyEqs     = Indices::numPrimaryEnergyVars};

    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state.
     *
     * If chemical equilibrium is not considered, we use the most
     * generic fluid state.
     */
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param priVars The primary variables
     * \param element The finite Element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     * \param temperatureIdx The temperature Index
     */
    Scalar getTemperature(const PrimaryVariables & priVars,
                          const Element & element,
                          const FVElementGeometry & fvGeometry,
                          const unsigned int scvIdx,
                          const Problem & problem,
                          const unsigned int temperatureIdx) const
    {
        // retrieve temperature from solution vector
        return priVars[temperature0Idx + temperatureIdx]; // could also be solid phase
    }

    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param priVars The primary Variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void updateTemperatures(FluidState & fluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars,
                            const Element & element,
                            const FVElementGeometry & fvGeometry,
                            const unsigned int scvIdx,
                            const Problem & problem)
    {
        assert(numPhases + 1 == numEnergyEqs);
        for(int phaseIdx=0; phaseIdx < numPhases; ++phaseIdx){
            // retrieve temperature from solution vector
            const Scalar T = priVars[temperature0Idx + phaseIdx];
            temperature_[phaseIdx]= T;
            fluidState.setTemperature(phaseIdx, T);
        }

        temperature_[sPhaseIdx] = priVars[temperature0Idx + sPhaseIdx];

      Valgrind::CheckDefined(temperature_);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     *
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void update(FluidState & fluidState,
                ParameterCache & paramCache,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx,
                const Problem & problem)
    {
        solidHeatCapacity_ =
            problem.spatialParams().solidHeatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidHeatCapacity_);

        for(int phaseIdx =0; phaseIdx<numPhases; ++phaseIdx){
            fluidThermalConductivity_[phaseIdx] =
          FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
        }
        Valgrind::CheckDefined(fluidThermalConductivity_);


        solidDensity_ =
                problem.spatialParams().solidDensity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidDensity_);

        solidThermalConductivity_ =
                problem.spatialParams().solidThermalConductivity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidThermalConductivity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            Valgrind::CheckDefined(h);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(kg K)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidHeatCapacity_; }

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperature(const unsigned int phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief Returns the total density of the given solid phase [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    { return solidDensity_; }

    /*!
     * \brief Returns the conductivity of the given solid phase [W/(m K)] in
     *        the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return solidThermalConductivity_; }

    /*!
     * \brief Returns the conductivity of the given fluid [W//m K)] in
     *        the sub-control volume.
     *
     *   \param phaseIdx The local index of the phases
     */
    Scalar thermalConductivity(const unsigned int phaseIdx) const
    {
        if(phaseIdx == wPhaseIdx or phaseIdx == nPhaseIdx )
            return fluidThermalConductivity_[phaseIdx];
        else if (phaseIdx == sPhaseIdx )
            return solidThermalConductivity_;
        else
            DUNE_THROW(Dune::NotImplemented,
                    "wrong index");
    }

    void checkDefinedTemp() const
    { Valgrind::CheckDefined(temperature_); }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(fluidThermalConductivity_);
        Valgrind::CheckDefined(solidThermalConductivity_);
        Valgrind::CheckDefined(solidDensity_);
        Valgrind::CheckDefined(solidHeatCapacity_);
    }

protected:
    Scalar temperature_[numPhases + 1];
    Scalar solidHeatCapacity_;
    Scalar solidDensity_;
    Scalar solidThermalConductivity_;
    Scalar fluidThermalConductivity_[numPhases];
};



/*!
 * \brief Contains the volume variables of the kinetic energy transfer
 *        module of the M-phase, N-component model.
 *        Specialization for the case of *2* energy balance equations.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/2>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { temperature0Idx = Indices::temperature0Idx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { sPhaseIdx = FluidSystem::sPhaseIdx };
    enum { numEnergyEqs     = Indices::numPrimaryEnergyVars};

    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state.
     *
     * If chemical equilibrium is not considered, we use the most
     * generic fluid state.
     */
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param priVars The primary variables
     * \param element The finite Element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     * \param temperatureIdx The temperature Index
     */
    Scalar getTemperature(const PrimaryVariables & priVars,
                          const Element & element,
                          const FVElementGeometry & fvGeometry,
                          const unsigned int scvIdx,
                          const Problem & problem,
                          const unsigned int temperatureIdx) const
    {
        // retrieve temperature from solution vector
        return priVars[temperature0Idx + temperatureIdx]; // could also be solid phase
    }

    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param priVars The primary Variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void updateTemperatures(FluidState & fluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars,
                            const Element & element,
                            const FVElementGeometry & fvGeometry,
                            const unsigned int scvIdx,
                            const Problem & problem)
    {
        assert(2 == numEnergyEqs);

        for(int energyEqIdx=0; energyEqIdx < numPhases; ++energyEqIdx){
            // retrieve temperature from solution vector
            const Scalar T = priVars[temperature0Idx + energyEqIdx];
            temperature_[energyEqIdx]= T;
        }

        fluidState.setTemperature(priVars[temperature0Idx]);

      Valgrind::CheckDefined(temperature_);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     *
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void update(FluidState & fluidState,
                ParameterCache & paramCache,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx,
                const Problem & problem)
    {
        solidHeatCapacity_ =
            problem.spatialParams().solidHeatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidHeatCapacity_);

        for(int phaseIdx =0; phaseIdx<numPhases; ++phaseIdx){
            fluidThermalConductivity_[phaseIdx] =
                    FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
        }
        Valgrind::CheckDefined(fluidThermalConductivity_);


        solidDensity_ =
                problem.spatialParams().solidDensity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidDensity_);

        solidThermalConductivity_ =
                problem.spatialParams().solidThermalConductivity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidThermalConductivity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            Valgrind::CheckDefined(h);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(kg K)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidHeatCapacity_; }

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperature(const unsigned int phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief Returns the total density of the given solid phase [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    { return solidDensity_; }

    /*!
     * \brief Returns the conductivity of the given solid phase [W/(m K)] in
     *        the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return solidThermalConductivity_; }

    /*!
     * \brief Returns the conductivity of the given fluid [W/(m K)] in
     *        the sub-control volume.
     *
     *   \param phaseIdx The local index of the phases
     */
    Scalar thermalConductivity(const unsigned int phaseIdx) const
    {
        if(phaseIdx == wPhaseIdx or phaseIdx == nPhaseIdx )
            return fluidThermalConductivity_[phaseIdx];
        else if (phaseIdx == sPhaseIdx )
            return solidThermalConductivity_;
        else
            DUNE_THROW(Dune::NotImplemented,
                    "wrong index");
    }

    void checkDefinedTemp() const
    { Valgrind::CheckDefined(temperature_); }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(fluidThermalConductivity_);
        Valgrind::CheckDefined(solidThermalConductivity_);
        Valgrind::CheckDefined(solidDensity_);
        Valgrind::CheckDefined(solidHeatCapacity_);
    }

protected:
    Scalar temperature_[numEnergyEqs];
    Scalar solidHeatCapacity_;
    Scalar solidDensity_;
    Scalar solidThermalConductivity_;
    Scalar fluidThermalConductivity_[numPhases];
};

} // end namespace

#endif
