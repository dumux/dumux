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
 * \ingroup
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_ENERGY_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/volumevariables.hh>

namespace Dumux {

// forward declaration
template <class IsothermalTraits, class Impl, bool enableEnergyBalance>
class EnergyVolumeVariablesImplementation;

/*!
 * \ingroup NIModel
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities. The volume variables base class
 *        is specialized for isothermal and non-isothermal models.
 */
template<class IsothermalTraits, class Impl>
using EnergyVolumeVariables = EnergyVolumeVariablesImplementation<IsothermalTraits,Impl, IsothermalTraits::ModelTraits::enableEnergyBalance()>;

/*!
 * \ingroup NIModel
 * \brief The isothermal base class
 */
template<class IsothermalTraits, class Impl>
class EnergyVolumeVariablesImplementation<IsothermalTraits, Impl, false>
{
    using Scalar = typename IsothermalTraits::PrimaryVariables::value_type;

public:
    using FluidState = typename IsothermalTraits::FluidState;
    using SolidState = typename IsothermalTraits::SolidState;
    using FluidSystem = typename IsothermalTraits::FluidSystem;

    //! The temperature is obtained from the problem as a constant for isothermal models
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateTemperature(const ElemSol& elemSol,
                           const Problem& problem,
                           const Element& element,
                           const Scv& scv,
                           FluidState& fluidState,
                           SolidState& solidState)
    {
        // retrieve temperature from solution vector, all phases have the same temperature
        Scalar T = problem.temperatureAtPos(scv.dofPosition());
        for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
        {
            fluidState.setTemperature(phaseIdx, T);
        }
        solidState.setTemperature(T);
    }

    template<class ElemSol, class Problem, class Element, class Scv>
    void updateSolidEnergyParams(const ElemSol &elemSol,
                                 const Problem& problem,
                                 const Element &element,
                                 const Scv &scv,
                                 SolidState & solidState)
    {}

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int phaseIdx)
    {
        return 0;
    }

};

//! The non-isothermal implicit volume variables base class
template<class IsothermalTraits, class Impl>
class EnergyVolumeVariablesImplementation<IsothermalTraits, Impl, true>
{
    using Scalar = typename IsothermalTraits::PrimaryVariables::value_type;
    using Idx = typename IsothermalTraits::ModelTraits::Indices;
    using ParentType = PorousMediumFlowVolumeVariables<IsothermalTraits>;

    static const int temperatureIdx = Idx::temperatureIdx;
    static const int numEnergyEq = IsothermalTraits::ModelTraits::numEnergyEq();

public:
    // export the fluidstate
    using FluidState = typename IsothermalTraits::FluidState;
    using SolidState = typename IsothermalTraits::SolidState;
    //! export the underlying fluid system
    using FluidSystem = typename IsothermalTraits::FluidSystem;
    using SolidSystem = typename IsothermalTraits::SolidSystem;

    //! The temperature is obtained from the problem as a constant for isothermal models
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateTemperature(const ElemSol& elemSol,
                           const Problem& problem,
                           const Element& element,
                           const Scv& scv,
                           FluidState& fluidState,
                           SolidState& solidState)
    {
        if (numEnergyEq == 1)
        {
            // retrieve temperature from solution vector, all phases have the same temperature
            const Scalar T = elemSol[scv.localDofIndex()][temperatureIdx];
            for(int phaseIdx=0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                fluidState.setTemperature(phaseIdx, T);
            }
            solidState.setTemperature(T);
        }

        else
        {
            for(int phaseIdx=0; phaseIdx < numEnergyEq-1; ++phaseIdx)
            {
                // retrieve temperatures from solution vector, phases might have different temperature
                const Scalar T = elemSol[scv.localDofIndex()][temperatureIdx + phaseIdx];
                fluidState.setTemperature(phaseIdx, T);
            }

            const Scalar solidTemperature = elemSol[scv.localDofIndex()][temperatureIdx+numEnergyEq-1];
            solidState.setTemperature(solidTemperature);
        }
    }

    template<class ElemSol, class Problem, class Element, class Scv>
    void updateSolidEnergyParams(const ElemSol &elemSol,
                                 const Problem& problem,
                                 const Element &element,
                                 const Scv &scv,
                                 SolidState & solidState)
    {

        Scalar solidHeatCapacity = SolidSystem::heatCapacity(solidState);
        solidState.setHeatCapacity(solidHeatCapacity);

        Scalar solidDensity = SolidSystem::density(solidState);
        solidState.setDensity(solidDensity);

        Scalar solidThermalConductivity = SolidSystem::thermalConductivity(solidState);
        solidState.setThermalConductivity(solidThermalConductivity);
    }

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(const int phaseIdx) const
    { return asImp_().fluidState().internalEnergy(phaseIdx); }

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(const int phaseIdx) const
    { return asImp_().fluidState().enthalpy(phaseIdx); }

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperatureSolid() const
    { return asImp_().solidState().temperature(); }


    /*!
     * \brief Returns the temperature of a fluid phase assuming thermal nonequilibrium
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperatureFluid(const int phaseIdx) const
    { return asImp_().fluidState().temperature(phaseIdx); }

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(kg K)]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return asImp_().solidState().heatCapacity(); }

    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    {  return  asImp_().solidState().density(); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of the solid phase in the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return asImp_().solidState().thermalConductivity(); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of a fluid phase in
     *        the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(asImp_().fluidState(), phaseIdx); }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

protected:
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }

};

} // end namespace Dumux

#endif
