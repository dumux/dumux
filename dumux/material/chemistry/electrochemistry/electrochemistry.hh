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
 *
 * \brief Electrochemical model for a fuel cell application.
 */
#ifndef ELECTRO_CHEM_HH
#define ELECTRO_CHEM_HH

#include <cmath>

#include <dumux/common/exceptions.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/fluidsystems/h2on2o2fluidsystem.hh>

namespace Dumux
{

/*!
 * \brief
 * The type of electrochemistry models implemented here (Ochs [2008] or Acosta et al. [2006])
 */
enum ElectroChemistryModel { Ochs, Acosta };

/*!
 * \brief
 * This class calculates source terms and current densities for fuel cells
 * with the electrochemical models suggested by Ochs [2008] or Acosta et al. [2006]
 */
template <class TypeTag, Dumux::ElectroChemistryModel electroChemistryModel>
class ElectroChemistry
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dumux::Constants<Scalar> Constant;

    typedef ElectroChemistry<TypeTag, electroChemistryModel> ThisType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        //indices of the phases
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum {
        //indices of the components
        wCompIdx = FluidSystem::wCompIdx, //major component of the liquid phase
        nCompIdx = FluidSystem::nCompIdx, //major component of the gas phase
        O2Idx = wCompIdx + 2
    };
    enum {
        //indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        switchIdx = Indices::switchIdx, //liquid saturation or mole fraction
        temperatureIdx = FluidSystem::numComponents //temperature
    };
    enum {
        //equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        contiH2OEqIdx = conti0EqIdx + wCompIdx,
        contiO2EqIdx = conti0EqIdx + wCompIdx + 2,
        energyEqIdx = FluidSystem::numComponents //energy equation
    };

public:
    /*!
    * \brief Calculates reaction sources with an electrochemical model approach.
    *
    * \param values The primary variable vector
    * \param volVars The volume variables
    *
    * For this method, the \a values parameter stores source values
    */
    static void reactionSource(PrimaryVariables &values,
                               const VolumeVariables &volVars)
    {
        static Scalar transportNumberH2O = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.TransportNumberH20);
        static Scalar maxIter = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.MaxIterations);
        static Scalar gridYMin = 0.0;
        static Scalar gridYMax = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.UpperRightY);
        static Scalar nCellsY = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.NumberOfCellsY);

        //initialise current density
        Scalar currentDensity = 0.0;

        //call internal method to calculate the current density
        currentDensity = calculateCurrentDensity_(volVars, maxIter);

        //correction to account for actually relevant reaction area
        //current density has to be devided by the half length of the box
        Scalar lengthBox = (gridYMax - gridYMin)/nCellsY;

        if(electroChemistryModel == ElectroChemistryModel::Acosta)
            currentDensity = currentDensity/lengthBox;
        else
            currentDensity = currentDensity*2/lengthBox;

        //conversion from [A/cm^2] to [A/m^2]
        currentDensity = currentDensity*10000;

        //calculation of flux terms with faraday equation
        values[contiH2OEqIdx] = currentDensity/(2*Constant::F);                  //reaction term in reaction layer
        values[contiH2OEqIdx] += currentDensity/Constant::F*transportNumberH2O;  //osmotic term in membrane
        values[contiO2EqIdx]  = -currentDensity/(4*Constant::F);                 //O2-equation
    }

protected:

    /*!
    * \brief Newton solver for calculation of the current density.
    */
    static Scalar calculateCurrentDensity_(const VolumeVariables &volVars, Scalar maxIter)
    {
        static Scalar specificResistance = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.SpecificResistance);
        static Scalar reversibleVoltage = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.ReversibleVoltage);
        static Scalar cellVoltage = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.CellVoltage);

        //initial guess for the current density and initial newton solver parameters
        Scalar currentDensity = reversibleVoltage - cellVoltage - 0.5;
        Scalar increment = 1e-4;
        Scalar deltaCurrentDensity = currentDensity*increment;
        Scalar deltaVoltage = 1.0;
        int iterations = 0;

        //Newton Solver for current Density
        while (std::abs(deltaVoltage) > 1e-6)
        {

            Scalar activationLosses        = calculateActivationLosses_(volVars, currentDensity);
            Scalar activationLossesNext    = calculateActivationLosses_(volVars, currentDensity+deltaCurrentDensity);
            Scalar concentrationLosses     = calculateConcentrationLosses_(volVars);
            Scalar activationLossesDiff    = activationLossesNext - activationLosses;
            Scalar sw                      = volVars.saturation(wPhaseIdx);

            if(electroChemistryModel == ElectroChemistryModel::Acosta)
            {
                // Acosta calculation
                deltaVoltage = currentDensity*specificResistance - reversibleVoltage + cellVoltage + activationLosses + concentrationLosses;

                currentDensity = currentDensity - (deltaVoltage*deltaCurrentDensity)/(deltaCurrentDensity*specificResistance + activationLossesDiff);

                activationLosses = calculateActivationLosses_(volVars, currentDensity);

                deltaVoltage = currentDensity*specificResistance - reversibleVoltage + cellVoltage + activationLosses + concentrationLosses;

                iterations++;

            }
            else
            {
                // Ochs & other calculation
                deltaVoltage = currentDensity*specificResistance/(1-sw) - reversibleVoltage + cellVoltage + activationLosses + concentrationLosses;

                currentDensity = currentDensity - (deltaVoltage*deltaCurrentDensity)/(deltaCurrentDensity*specificResistance/(1-sw) + activationLossesDiff);

                activationLosses = calculateActivationLosses_(volVars, currentDensity);

                deltaVoltage = currentDensity*specificResistance/(1-sw) - reversibleVoltage + cellVoltage + activationLosses + concentrationLosses;

                iterations++;
            }

            if(iterations >= maxIter)
            {
                DUNE_THROW(Dumux::NumericalProblem, "Newton solver for electrochemistry didn't converge");
            }
        }

        return currentDensity;
    }

private:

    /*!
    * \brief Calculation of the activation losses
    */
    static Scalar calculateActivationLosses_(const VolumeVariables &volVars, const Scalar currentDensity)
    {
        static Scalar refO2PartialPressure = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.RefO2PartialPressure);
        static Scalar numElectrons = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.NumElectrons);
        static Scalar transferCoefficient = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.TransferCoefficient);

        //Saturation sw for Acosta calculation
        Scalar sw = volVars.saturation(wPhaseIdx);
        //Calculate prefactor
        Scalar preFactor = Constant::R*volVars.fluidState().temperature()/transferCoefficient/Constant::F/numElectrons;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(nPhaseIdx) * volVars.fluidState().moleFraction(nPhaseIdx, O2Idx);

        Scalar losses = 0.0;
        //Calculate activation losses
        if(electroChemistryModel == ElectroChemistryModel::Acosta)
        {
            losses = preFactor
                            *(  std::log(std::abs(currentDensity)/std::abs(exchangeCurrentDensity_(volVars)))
                            - std::log(pO2/refO2PartialPressure)
                            - std::log(1 - sw)
                            );
        }
        else
        {
            losses = preFactor
            *(  std::log(std::abs(currentDensity)/std::abs(exchangeCurrentDensity_(volVars)))
                - std::log(pO2/refO2PartialPressure)
                );
        }
        return losses;
    }


    /*!
    * \brief Calculation of concentration losses.
    */
    static Scalar calculateConcentrationLosses_(const VolumeVariables &volVars)
    {
        static Scalar pO2Inlet = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.pO2Inlet);
        static Scalar numElectrons = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.NumElectrons);
        static Scalar transferCoefficient = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.TransferCoefficient);

        //Calculate preFactor
        Scalar preFactor = Constant::R*volVars.temperature()/transferCoefficient/Constant::F/numElectrons;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(nPhaseIdx) * volVars.fluidState().moleFraction(nPhaseIdx, O2Idx);

        Scalar losses = 0.0;
        //Calculate concentration losses
        if(electroChemistryModel == ElectroChemistryModel::Acosta)
        {
            losses = -1.0*preFactor*(transferCoefficient/2)*std::log(pO2/pO2Inlet);
        }else
        {
            // +1 is the Nernst part of the equation
            losses = -1.0*preFactor*(transferCoefficient/2+1)*std::log(pO2/pO2Inlet);
        }

        return losses;
    }


    /*!
    * \brief Calculation of the exchange current density.
    */
    static Scalar exchangeCurrentDensity_(const VolumeVariables &volVars)
    {
        static Scalar activationBarrier = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.ActivationBarrier);
        static Scalar surfaceIncreasingFactor = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.SurfaceIncreasingFactor);
        static Scalar refTemperature = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.RefTemperature);
        static Scalar refCurrentDensity = GET_RUNTIME_PARAM(TypeTag, Scalar, ElectroChemistry.RefCurrentDensity);

        Scalar T = volVars.fluidState().temperature();
        Scalar refExchangeCurrentDensity = -1.0
                            * refCurrentDensity
                            * surfaceIncreasingFactor
                            * std::exp(-1.0 * activationBarrier / Constant::R * (1/T-1/refTemperature));

        return refExchangeCurrentDensity;
    }
};

}// end namespace

#endif
