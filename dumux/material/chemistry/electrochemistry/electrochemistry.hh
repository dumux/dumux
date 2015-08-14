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

#include <dumux/common/exceptions.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/fluidsystems/h2on2o2fluidsystem.hh>

#include <cmath>
#include <iostream>
namespace Dumux
{
    
static const int Ochs = 1;
static const int Acosta = 2;

    /*!
     * \brief 
     * This class calculates source terms and current densities for fuel cells 
     * with the electrochemical models suggested by Ochs [2008] or Acosta [2006]
     */
    template <class TypeTag, const int electroChemApproach>
    class ElectroChemestry
    {
    protected:
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
        typedef typename GridView::template Codim<0>::Entity Element;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
        
        typedef Dumux::Constants<Scalar> Constant;

        typedef ElectroChemestry<TypeTag,electroChemApproach> ThisType;
        
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        
        enum { //indices of the phases
                wPhaseIdx = Indices::wPhaseIdx,
                nPhaseIdx = Indices::nPhaseIdx,
                
                //indices of the components
                wCompIdx = FluidSystem::wCompIdx, //major component of the liquid phase
                nCompIdx = FluidSystem::nCompIdx, //major component of the gas phase
                O2Idx = wCompIdx + 2,

                //indices of the primary variables
                pressureIdx = Indices::pressureIdx, //gas-phase pressure
                switchIdx = Indices::switchIdx, //liquid saturation or mole fraction
                temperatureIdx = FluidSystem::numComponents, //temperature
                
                //equation indices
                conti0EqIdx = Indices::conti0EqIdx,
                contiH2OEqIdx = conti0EqIdx + wCompIdx,
                contiO2EqIdx = conti0EqIdx + wCompIdx + 2,
                energyEqIdx = FluidSystem::numComponents, //energy equation
            };
        
    public:
        ElectroChemestry()
        {
            gridYMin_ = 0.0;
            gridYMax_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.UpperRightY);
            nCellsY_  = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.NumberOfCellsY);
            
            maxIter_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.MaxIterations);
            
            cellVoltage_        	 = 	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.CellVoltage);
            thermoneutralVoltage_ 	 = 	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.ThermoneutralVoltage);
            reversibleVoltage_  	 =	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.ReversibleVoltage);
            specificResistance_ 	 = 	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.SpecificResistance);
            
            transferCoefficient_	 =  GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.TransferCoefficient);
            numElectrons_			 =  GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.NumElectrons);
            refO2PartialPressure_	 =  GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.RefO2PartialPressure);
            transportNumberH2O_		 =  GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.TransportNumberH20);
            pO2Inlet_				 =  GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.pO2Inlet); 
            refCurrentDensity_ 		 = 	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.RefCurrentDensity);
            refTemperature_			 =	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.RefTemperature);
            activationBarrier_   	 =	GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.ActivationBarrier);
            surfaceIncreasingFactor_  = GET_RUNTIME_PARAM(TypeTag, Scalar, FuelCell.SurfaceIncreasingFactor);
            eps_ = 1e-6;
        }
                    
        /*!
        * \brief Calculates reaction sources with an electrochemical model approach.
        *
        * \param q The primary variable vector
        * \param volVars The volume variables
        *
        * For this method, the \a q parameter stores primary
        * variables.
        */
        void reactionSource(PrimaryVariables &q,
                            const VolumeVariables &volVars)
        {
            //initialise current density
            Scalar currentDensity = 0.0;
            
            //call internal method to calculate the current density
            currentDensity = calculateCurrentDensity_(volVars, maxIter_);
                        
            //correction to account for actually relevant reaction area
            //current density has to be devided by the half length of the box
            //see diploma thesis Lena Walter
            Scalar lengthBox= (gridYMax_ - gridYMin_)/nCellsY_;
            
            if(electroChemApproach==2)
                currentDensity = currentDensity/lengthBox; //ACOSTA
            else
                currentDensity = currentDensity*2/lengthBox;
            
            //conversion from [A/cm^2] to [A/m^2]
            currentDensity = currentDensity*10000; 
            
            //calculation of flux terms with faraday equation
            q[contiH2OEqIdx] = currentDensity/(2*Constant::F);                  //reaction term in reaction layer
            q[contiH2OEqIdx] += currentDensity/Constant::F*transportNumberH2O_; //osmotic term in membrane
            q[contiO2EqIdx]  = -currentDensity/(4*Constant::F);                 //O2-equation
        }
        
protected:
    
    /*!
    * \brief Newton solver for calculation of the current density.
    */
    Scalar calculateCurrentDensity_(const VolumeVariables &volVars, Scalar maxIter_)
    {
        //initial guess for the current density and initial newton solver parameters
        Scalar currentDensity = reversibleVoltage_ - cellVoltage_ - 0.5; 
        Scalar increment = 1e-4;
        Scalar deltaCurrentDensity = currentDensity*increment;
        Scalar deltaVoltage = 1.0;
        int iterations = 0;
        
        //Newton Solver for current Density
        while (absolute_(deltaVoltage)>eps_) 
        {
            
            Scalar activationLosses 		= 	calculateActivationLosses_(volVars, currentDensity);
            Scalar activationLossesNext 	= 	calculateActivationLosses_(volVars, currentDensity+deltaCurrentDensity);
            Scalar concentrationLosses 		= 	calculateConcentrationLosses_(volVars);
            Scalar activationLossesDiff 	= 	activationLossesNext - activationLosses;
            Scalar sw						= 	volVars.saturation(wPhaseIdx);
                                
            if(electroChemApproach==2)
            {
                //Acosta calculation
                deltaVoltage = currentDensity*specificResistance_ - reversibleVoltage_ + cellVoltage_ + activationLosses + concentrationLosses;
                                
                currentDensity = currentDensity - (deltaVoltage*deltaCurrentDensity)/(deltaCurrentDensity*specificResistance_ + activationLossesDiff);
                                                                
                activationLosses = calculateActivationLosses_(volVars, currentDensity);
                
                deltaVoltage = currentDensity*specificResistance_ - reversibleVoltage_ + cellVoltage_ + activationLosses + concentrationLosses;
                                
                iterations++;
            
            }
            else
            {
                //Ochs calculation
                deltaVoltage = currentDensity*specificResistance_/(1-sw) - reversibleVoltage_ + cellVoltage_ + activationLosses + concentrationLosses;

                currentDensity = currentDensity - (deltaVoltage*deltaCurrentDensity)/(deltaCurrentDensity*specificResistance_/(1-sw) + activationLossesDiff);
                
                activationLosses = calculateActivationLosses_(volVars, currentDensity);

                deltaVoltage = currentDensity*specificResistance_/(1-sw) - reversibleVoltage_ + cellVoltage_ + activationLosses + concentrationLosses;

                iterations++;
            }
            
            if(iterations >= maxIter_)
            {
                DUNE_THROW(Dumux::NumericalProblem, "Newton solver for electrochemistry didn't converge");
            }
        }
        
        return currentDensity;
    }
    
    /*!
    * \brief Calculation of the activation losses
    */		
    Scalar calculateActivationLosses_(const VolumeVariables &volVars, const Scalar currentDensity)
    {
        //Saturation sw for Acosta calculation
        Scalar sw = volVars.saturation(wPhaseIdx); 
        //Calculate prefactor
        Scalar preFactor = Constant::R*volVars.fluidState().temperature()/transferCoefficient_/Constant::F/numElectrons_;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(nPhaseIdx) * volVars.fluidState().moleFraction(nPhaseIdx, O2Idx); 			
        
        Scalar losses = 0.0;
        //Calculate activation losses
        if(electroChemApproach==2) 
        {
            losses = preFactor
                            *(  log(absolute_(currentDensity)/absolute_(exchangeCurrentDensity_(volVars)))
                            - log(pO2/refO2PartialPressure_)
                            - log(1 - sw)//This Term only for Acosta calculation
                            );
        } 
        else 
        {
            losses = preFactor
            *(  log(absolute_(currentDensity)/absolute_(exchangeCurrentDensity_(volVars)))
                - log(pO2/refO2PartialPressure_)
                );
        }
        return losses;
    }
    
    
    /*!
    * \brief Calculation of concentration losses.
    */
    Scalar calculateConcentrationLosses_(const VolumeVariables &volVars)
    {
        //Calculate preFactor
        Scalar preFactor = Constant::R*volVars.temperature()/transferCoefficient_/Constant::F/numElectrons_;
        //Get partial pressure of O2 in the gas phase
        Scalar pO2 = volVars.pressure(nPhaseIdx) * volVars.fluidState().moleFraction(nPhaseIdx, O2Idx); 
        
        Scalar losses = 0.0;
        //Calculate concentration losses
        if(electroChemApproach==2) 
        {
            //Acosta
            losses = -1.0*preFactor*(transferCoefficient_/2)*log(pO2/pO2Inlet_);
        }else
        {
            //Ochs (+1 is the Nernst part of the equation)
            losses = -1.0*preFactor*(transferCoefficient_/2+1)*log(pO2/pO2Inlet_);
        }

        return losses;
    }
    
        
    /*!
    * \brief Calculation of the exchange current density.
    */
    Scalar exchangeCurrentDensity_(const VolumeVariables &volVars)
    {
        Scalar T = volVars.fluidState().temperature();
        Scalar refExchangeCurrentDensity = -1.0 
                            * refCurrentDensity_
                            * surfaceIncreasingFactor_
                            * exp(-1.0* activationBarrier_ /Constant::R*(1/T-1/refTemperature_));
        
        return refExchangeCurrentDensity;
    }
    
    /*!
    * \brief Internal method to return the absolute value of a scalar.
    */
    Scalar absolute_(Scalar x)
    {
        if(x<0.0)
        {
            return x*(-1);
        }
        else return x;
    }


    Scalar  eps_; 
    Scalar	gridYMin_;
    Scalar	gridYMax_;
    Scalar	nCellsY_;
    Scalar	maxIter_;
    Scalar	cellVoltage_;
    Scalar	thermoneutralVoltage_;
    Scalar	reversibleVoltage_;
    Scalar	specificResistance_;
    Scalar	transferCoefficient_;
    Scalar	numElectrons_;
    Scalar	refO2PartialPressure_;
    Scalar	transportNumberH2O_;
    Scalar	pO2Inlet_; 
    Scalar	refCurrentDensity_;
    Scalar	refTemperature_;
    Scalar	activationBarrier_;
    Scalar	surfaceIncreasingFactor_;
}; 
}// end namespace
#endif
	

