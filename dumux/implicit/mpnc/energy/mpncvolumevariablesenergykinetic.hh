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

#include <dumux/implicit/mpnc/mpncvolumevariables.hh>
#include "mpncvtkwriterenergykinetic.hh"

namespace Dumux
{
/*!
 * \brief Contains the volume variables of the kinetic energy transfer
 *        module of the M-phase, N-component model.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { temperature0Idx = Indices::temperature0Idx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { sPhaseIdx = FluidSystem::sPhaseIdx };
    enum { numEnergyEqs     = Indices::NumPrimaryEnergyVars};
    enum { dim = GridView::dimension};

    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state.
     *
     * If chemical equilibrium is not considered, we use the most
     * generic fluid state.
     */
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param The primary variables
     * \param element The finit Element
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
        heatCapacity_ =
            problem.spatialParams().heatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);

        for(int phaseIdx =0; phaseIdx<numPhases; ++phaseIdx){
            fluidThermalConductivity_[phaseIdx] =
	      FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
        }
        Valgrind::CheckDefined(fluidThermalConductivity_);


        soilDensity_ =
                problem.spatialParams().soilDensity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(soilDensity_);

        soilThermalConductivity_ =
                problem.spatialParams().soilThermalConductivity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(soilThermalConductivity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            Valgrind::CheckDefined(h);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(K m^3)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; }

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperature(const unsigned int phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief Returns the total density of the given soil [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar soilDensity() const
    { return soilDensity_; }

    /*!
     * \brief Returns the conductivity of the given soil [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar soilThermalConductivity() const
    { return soilThermalConductivity_; }

    /*!
     * \brief Returns the conductivity of the given fluid [kg / m^3] in
     *        the sub-control volume.
     *
     *   \param phaseIdx The local index of the phases
     */
    Scalar thermalConductivity(const unsigned int phaseIdx) const
    {
        if(phaseIdx == wPhaseIdx or phaseIdx == nPhaseIdx )
            return fluidThermalConductivity_[phaseIdx];
        else if (phaseIdx == sPhaseIdx )
            return soilThermalConductivity_;
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
        Valgrind::CheckDefined(soilThermalConductivity_);
        Valgrind::CheckDefined(soilDensity_);
        Valgrind::CheckDefined(heatCapacity_);
    }

protected:
    Scalar temperature_[numPhases + 1];
    Scalar heatCapacity_;
    Scalar soilDensity_;
    Scalar soilThermalConductivity_;
    Scalar fluidThermalConductivity_[numPhases];
};

} // end namepace

#endif
