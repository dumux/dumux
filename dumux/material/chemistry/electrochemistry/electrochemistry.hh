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
 * \ingroup Chemistry
 * \brief Electrochemical model for a fuel cell application.
 */
#ifndef DUMUX_ELECTROCHEMISTRY_HH
#define DUMUX_ELECTROCHEMISTRY_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/discretization/method.hh>
#include <dumux/material/constants.hh>

#include <dumux/material/fluidsystems/h2on2o2.hh>

namespace Dumux {

/*!
 * \ingroup Chemistry
 * \brief The type of electrochemistry models implemented here (Ochs (2008) \cite ochs2008 or Acosta et al. (2006) \cite A3:acosta:2006 )
 */
enum ElectroChemistryModel { Ochs, Acosta };

/*!
 * \ingroup Chemistry
 * \brief This class calculates source terms and current densities for fuel cells
 * with the electrochemical models suggested by Ochs (2008) \cite ochs2008 or Acosta et al. (2006) \cite A3:acosta:2006
 * \todo TODO: Scalar type should be extracted from VolumeVariables!
 * \todo TODO: This shouldn't depend on grid and discretization!!
 */
template <class Scalar, class Indices, class FluidSystem, class GridGeometry, ElectroChemistryModel electroChemistryModel>
class ElectroChemistry
{

    using Constant = Dumux::Constants<Scalar>;

    enum
    {
        //indices of the primary variables
        pressureIdx = Indices::pressureIdx, // gas-phase pressure
        switchIdx = Indices::switchIdx,     // liquid saturation or mole fraction
    };
    enum
    {
        //equation indices
        contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiO2EqIdx = Indices::conti0EqIdx + FluidSystem::O2Idx,
    };

    enum
    {
        // phase indices
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx
    };

    using GridView = typename GridGeometry::GridView;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    using GlobalPosition = typename Dune::FieldVector<typename GridView::ctype, GridView::dimensionworld>;

public:
    /*!
     * \brief Calculates reaction sources with an electrochemical model approach.
     *
     * \param values The primary variable vector
     * \param currentDensity The current density
     * \param paramGroup The group containing the needed parameters
     *
     * For this method, the \a values parameter stores source values
     */
    template<class SourceValues>
    static void reactionSource(SourceValues &values, Scalar currentDensity,
                               const std::string& paramGroup = "")
    {
        // correction to account for actually relevant reaction area
        // current density has to be divided by the half length of the box
        // \todo Do we have multiply with the electrochemically active surface area (ECSA) here instead?
        static Scalar gridYMax = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight")[1];
        static Scalar nCellsY = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.Cells")[1];

        // Warning: This assumes the reaction layer is always just one cell (cell-centered) or half a box (box) thick
        const auto lengthBox = gridYMax/nCellsY;
        if (isBox)
            currentDensity *= 2.0/lengthBox;
        else
            currentDensity *= 1.0/lengthBox;

        static Scalar transportNumberH2O = getParam<Scalar>("ElectroChemistry.TransportNumberH20");

        //calculation of flux terms with Faraday equation
        values[contiH2OEqIdx] = currentDensity/(2*Constant::F);                  //reaction term in reaction layer
        values[contiH2OEqIdx] += currentDensity/Constant::F*transportNumberH2O;  //osmotic term in membrane
        values[contiO2EqIdx]  = -currentDensity/(4*Constant::F);                 //O2-equation
    }

    /*!
     * \brief Newton solver for calculation of the current density.
     * \param volVars The volume variables
     * \returns The current density in A/m^2
     */
    template<class VolumeVariables>
    static Scalar calculateCurrentDensity(const VolumeVariables &volVars)
    {

        static int maxIter = getParam<int>("ElectroChemistry.MaxIterations");
        static Scalar specificResistance = getParam<Scalar>("ElectroChemistry.SpecificResistance");
        static Scalar reversibleVoltage = getParam<Scalar>("ElectroChemistry.ReversibleVoltage");
        static Scalar cellVoltage = getParam<Scalar>("ElectroChemistry.CellVoltage");

        //initial guess for the current density and initial newton solver parameters
        Scalar currentDensity = reversibleVoltage - cellVoltage - 0.5;
        Scalar increment = 1e-4;
        Scalar deltaCurrentDensity = currentDensity*increment;
        Scalar deltaVoltage = 1.0;
        int iterations = 0;

        //Newton Solver for current Density
        using std::abs;
        while (abs(deltaVoltage) > 1e-6)
        {

            Scalar activationLosses        = calculateActivationLosses_(volVars, currentDensity);
            Scalar activationLossesNext    = calculateActivationLosses_(volVars, currentDensity+deltaCurrentDensity);
            Scalar concentrationLosses     = calculateConcentrationLosses_(volVars);
            Scalar activationLossesDiff    = activationLossesNext - activationLosses;
            Scalar sw                      = volVars.saturation(liquidPhaseIdx);

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
                DUNE_THROW(NumericalProblem, "Newton solver for electrochemistry didn't converge");
            }
        }

        //conversion from [A/cm^2] to [A/m^2]
        return currentDensity*10000;
    }

private:

    /*!
     * \brief Calculation of the activation losses
     * \param volVars The volume variables
     * \param currentDensity The current density
     */
    template<class VolumeVariables>
    static Scalar calculateActivationLosses_(const VolumeVariables &volVars, const Scalar currentDensity)
    {
        static Scalar refO2PartialPressure = getParam<Scalar>("ElectroChemistry.RefO2PartialPressure");
        static Scalar numElectrons = getParam<Scalar>("ElectroChemistry.NumElectrons");
        static Scalar transferCoefficient = getParam<Scalar>("ElectroChemistry.TransferCoefficient");

        //Saturation sw for Acosta calculation
        Scalar sw = volVars.saturation(liquidPhaseIdx);
        //Calculate prefactor
        Scalar preFactor = Constant::R*volVars.fluidState().temperature()/transferCoefficient/Constant::F/numElectrons;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(gasPhaseIdx) * volVars.fluidState().moleFraction(gasPhaseIdx, FluidSystem::O2Idx);

        Scalar losses = 0.0;
        //Calculate activation losses
        using std::log;
        using std::abs;
        if(electroChemistryModel == ElectroChemistryModel::Acosta)
        {
            losses = preFactor
                            *(  log(abs(currentDensity)/abs(exchangeCurrentDensity_(volVars)))
                            - log(pO2/refO2PartialPressure)
                            - log(1 - sw)
                            );
        }
        else
        {
            losses = preFactor
            *(  log(abs(currentDensity)/abs(exchangeCurrentDensity_(volVars)))
                - log(pO2/refO2PartialPressure)
                );
        }
        return losses;
    }


    /*!
     * \brief Calculation of concentration losses.
     * \param volVars The volume variables
     */
    template<class VolumeVariables>
    static Scalar calculateConcentrationLosses_(const VolumeVariables &volVars)
    {
        static Scalar pO2Inlet = getParam<Scalar>("ElectroChemistry.pO2Inlet");
        static Scalar numElectrons = getParam<Scalar>("ElectroChemistry.NumElectrons");
        static Scalar transferCoefficient = getParam<Scalar>("ElectroChemistry.TransferCoefficient");

        //Calculate preFactor
        Scalar preFactor = Constant::R*volVars.temperature()/transferCoefficient/Constant::F/numElectrons;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(gasPhaseIdx) * volVars.fluidState().moleFraction(gasPhaseIdx, FluidSystem::O2Idx);

        Scalar losses = 0.0;
        //Calculate concentration losses
        using std::log;
        if(electroChemistryModel == ElectroChemistryModel::Acosta)
        {
            losses = -1.0*preFactor*(transferCoefficient/2)*log(pO2/pO2Inlet);
        }else
        {
            // +1 is the Nernst part of the equation
            losses = -1.0*preFactor*(transferCoefficient/2+1)*log(pO2/pO2Inlet);
        }

        return losses;
    }


    /*!
     * \brief Calculation of the exchange current density.
     * \param volVars The volume variables
     */
    template<class VolumeVariables>
    static Scalar exchangeCurrentDensity_(const VolumeVariables &volVars)
    {
        using std::exp;
        static Scalar activationBarrier = getParam<Scalar>("ElectroChemistry.ActivationBarrier");
        static Scalar surfaceIncreasingFactor = getParam<Scalar>("ElectroChemistry.SurfaceIncreasingFactor");
        static Scalar refTemperature = getParam<Scalar>("ElectroChemistry.RefTemperature");
        static Scalar refCurrentDensity = getParam<Scalar>("ElectroChemistry.RefCurrentDensity");

        Scalar T = volVars.fluidState().temperature();
        Scalar refExchangeCurrentDensity = -1.0
                            * refCurrentDensity
                            * surfaceIncreasingFactor
                            * exp(-1.0 * activationBarrier / Constant::R * (1/T-1/refTemperature));

        return refExchangeCurrentDensity;
    }
};

} // end namespace Dumux

#endif
