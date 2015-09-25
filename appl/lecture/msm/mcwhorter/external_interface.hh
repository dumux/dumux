// $Id: pcm_parameters.hh 1329 $
/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef PCM_PARAMETERS_HH
#define PCM_PARAMETERS_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>

/** \todo Please doc me! */
namespace Dumux
{

class InterfaceSoilProperties
//Class for Interface Soil Properties
//Integrate Parameters for Probabilistic Collocation Method (iPCM).
//Contained design and uncertainty parameters
{
public:
    //Interface Soil Properties (ISP):
    float ISP_Permeability; // Permeability
    float ISP_FinePermeability; // Fine Permeability
    float ISP_CoarsePermeability; // Coarse Permeability
    float ISP_Porosity; // Porosity
    float ISP_FinePorosity; // Fine Porosity
    float ISP_CoarsePorosity; // Coarse Porosity
    float ISP_LeakageWellPermeability; // Leakage Well Permeability
    float ISP_LongitudinalDispersivity; //Longitudinal dispersivity
    float ISP_TransverseDispersivity; //Transverse dispersivity
    float ISP_BrooksCoreyLambda;
    float ISP_BrooksCoreyEntryPressure;
    float ISP_FineBrooksCoreyLambda;
    float ISP_FineBrooksCoreyEntryPressure;
    float ISP_CoarseBrooksCoreyLambda;
    float ISP_CoarseBrooksCoreyEntryPressure;
    float ISP_MinCapillaryPressure;
    float ISP_MaxCapillaryPressure;
    float ISP_ResidualSaturationWetting;
    float ISP_ResidualSaturationNonWetting;
    float ISP_FineResidualSaturationWetting;
    float ISP_FineResidualSaturationNonWetting;
    float ISP_CoarseResidualSaturationWetting;
    float ISP_CoarseResidualSaturationNonWetting;

    InterfaceSoilProperties(const char* isp_filename)
    //Initialization of ISP Parameters
    {
        using namespace std;
        std::cout
                << "-----> ISP: Interface Soil Properties Initialization ...\n";
        //ISP input file defenition
        ifstream input;

        //ISP file check
        input.open(isp_filename);
        if (!input)
        {
            cout << "\n";
            cout << "-----> ISP: Fatal error! - Data read \n";
            cout << "-----> ISP: Could not open the input data file: \""
                    << isp_filename << "\n";
        }

        //iPCM input file reading
        char reader[100]; // variable for input value
        while (!input.eof())
        {
            input >> reader;
            //if (reader==string("<SoilProperties>"))
            //cout << "-----> ISP: Soil Properties reading ... \n";
            //ISP perameters initialization:
            if (reader == string("<Permeability>"))
            {
                input >> reader;
                ISP_Permeability = atof(reader);
                cout << "-----> ISP: Permeability: " << ISP_Permeability
                        << "\n";
            }
            if (reader == string("<FinePermeability>"))
            {
                input >> reader;
                ISP_FinePermeability = atof(reader);
                cout << "-----> ISP: Fine permeability: "
                        << ISP_FinePermeability << "\n";
            }
            if (reader == string("<CoarsePermeability>"))
            {
                input >> reader;
                ISP_CoarsePermeability = atof(reader);
                cout << "-----> ISP: Coarse permeability: "
                        << ISP_CoarsePermeability << "\n";
            }
            if (reader == string("<Porosity>"))
            {
                input >> reader;
                ISP_Porosity = atof(reader);
                cout << "-----> ISP: Porosity: " << ISP_Porosity << "\n";
            }
            if (reader == string("<FinePorosity>"))
            {
                input >> reader;
                ISP_FinePorosity = atof(reader);
                cout << "-----> ISP: Fine porosity: " << ISP_FinePorosity
                        << "\n";
            }
            if (reader == string("<CoarsePorosity>"))
            {
                input >> reader;
                ISP_CoarsePorosity = atof(reader);
                cout << "-----> ISP: Coarse porosity: " << ISP_CoarsePorosity
                        << "\n";
            }
            if (reader == string("<LeakageWellPermeability>"))
            {
                input >> reader;
                ISP_LeakageWellPermeability = atof(reader);
                cout << "-----> ISP: Leakage Well Permeability: "
                        << ISP_LeakageWellPermeability << "\n";
            }
            if (reader == string("<LongitudinalDispersivity>"))
            {
                input >> reader;
                ISP_LongitudinalDispersivity = atof(reader);
                cout << "-----> ISP: Longitudinal dispersivity: "
                        << ISP_LongitudinalDispersivity << "\n";
            }
            if (reader == string("<TransverseDispersivity>"))
            {
                input >> reader;
                ISP_TransverseDispersivity = atof(reader);
                cout << "-----> ISP: Transverse dispersivity: "
                        << ISP_TransverseDispersivity << "\n";
            }
            if (reader == string("<BrooksCoreyLambda>"))
            {
                input >> reader;
                ISP_BrooksCoreyLambda = atof(reader);
                cout << "-----> ISP: Brooks-Corey lambda: "
                        << ISP_BrooksCoreyLambda << "\n";
            }
            if (reader == string("<BrooksCoreyEntryPressure>"))
            {
                input >> reader;
                ISP_BrooksCoreyEntryPressure = atof(reader);
                cout << "-----> ISP: Brooks-Corey entry pressure: "
                        << ISP_BrooksCoreyEntryPressure << "\n";
            }

            if (reader == string("<FineBrooksCoreyLambda>"))
            {
                input >> reader;
                ISP_FineBrooksCoreyLambda = atof(reader);
                cout << "-----> ISP: Brooks-Corey lambda, fine: "
                        << ISP_FineBrooksCoreyLambda << "\n";
            }
            if (reader == string("<FineBrooksCoreyEntryPressure>"))
            {
                input >> reader;
                ISP_FineBrooksCoreyEntryPressure = atof(reader);
                cout << "-----> ISP: Brooks-Corey entry pressure, fine: "
                        << ISP_FineBrooksCoreyEntryPressure << "\n";
            }
            if (reader == string("<CoarseBrooksCoreyLambda>"))
            {
                input >> reader;
                ISP_CoarseBrooksCoreyLambda = atof(reader);
                cout << "-----> ISP: Brooks-Corey lambda, coarse: "
                        << ISP_CoarseBrooksCoreyLambda << "\n";
            }
            if (reader == string("<CoarseBrooksCoreyEntryPressure>"))
            {
                input >> reader;
                ISP_CoarseBrooksCoreyEntryPressure = atof(reader);
                cout << "-----> ISP: Brooks-Corey entry pressure, coarse: "
                        << ISP_CoarseBrooksCoreyEntryPressure << "\n";
            }
            if (reader == string("<MinCapillaryPressure>"))
            {
                input >> reader;
                ISP_MinCapillaryPressure = atof(reader);
                cout << "-----> ISP: Minimal capillary pressure: "
                        << ISP_MinCapillaryPressure << "\n";
            }
            if (reader == string("<MaxCapillaryPressure>"))
            {
                input >> reader;
                ISP_MaxCapillaryPressure = atof(reader);
                cout << "-----> ISP: Maximal capillary pressure: "
                        << ISP_MaxCapillaryPressure << "\n";
            }
            if (reader == string("<ResidualSaturationWetting>"))
            {
                input >> reader;
                ISP_ResidualSaturationWetting = atof(reader);
                cout << "-----> ISP: Residual saturation wetting phase: "
                        << ISP_ResidualSaturationWetting << "\n";
            }
            if (reader == string("<ResidualSaturationNonWetting>"))
            {
                input >> reader;
                ISP_ResidualSaturationNonWetting = atof(reader);
                cout << "-----> ISP: Residual saturation nonwetting phase: "
                        << ISP_ResidualSaturationNonWetting << "\n";
            }
            if (reader == string("<FineResidualSaturationWetting>"))
            {
                input >> reader;
                ISP_FineResidualSaturationWetting = atof(reader);
                cout << "-----> ISP: Residual saturation wetting phase, fine: "
                        << ISP_FineResidualSaturationWetting << "\n";
            }
            if (reader == string("<FineResidualSaturationNonWetting>"))
            {
                input >> reader;
                ISP_FineResidualSaturationNonWetting = atof(reader);
                cout << "-----> ISP: Residual saturation nonwetting phase, fine: "
                        << ISP_FineResidualSaturationNonWetting << "\n";
            }
            if (reader == string("<CoarseResidualSaturationWetting>"))
            {
                input >> reader;
                ISP_CoarseResidualSaturationWetting = atof(reader);
                cout << "-----> ISP: Residual saturation wetting phase, coarse: "
                        << ISP_CoarseResidualSaturationWetting << "\n";
            }
            if (reader == string("<CoarseResidualSaturationNonWetting>"))
            {
                input >> reader;
                ISP_CoarseResidualSaturationNonWetting = atof(reader);
                cout << "-----> ISP: Residual saturation nonwetting phase, coarse: "
                        << ISP_CoarseResidualSaturationNonWetting << "\n";
            }
        }
        input.close();
    }

};

class InterfaceFluidProperties
//Class for Interface Fluid Properties
//Integrate Parameters for Probabilistic Collocation Method (iPCM).
//Contained design and uncertainty parameters
{
public:
    //Interface Fluid Properties (IFP):
    float IFP_GasDiffCoeff;// Gas Diffusion Coefficient
    float IFP_CO2ResidSat;// Residual Saturation of CO2
    float IFP_MolecularDiffusionCoefficient;
    float IFP_ViscosityWettingFluid;
    float IFP_ViscosityNonWettingFluid;

    InterfaceFluidProperties(const char* ifp_filename)
    //Initialization of IFP Parameters
    {
        using namespace std;
        std::cout
                << "-----> IFP: Interface Fluid Properties Initialization ...\n";
        //IFP input file defenition
        ifstream input;

        //IFP file check
        input.open(ifp_filename);
        if (!input)
        {
            cout << "\n";
            cout << "-----> IFP: Fatal error! - Data read \n";
            cout << "-----> IFP: Could not open the input data file: \""
                    << ifp_filename << "\n";
        }

        //iPCM input file reading
        char reader[100]; // variable for input value
        //double K;
        while (!input.eof())
        {
            input >> reader;
            //if (reader==string("<FluidProperties>"))
            //cout << "-----> IFP: Fluid Properties reading ... \n";
            //IFP perameters initialization:
            if (reader == string("<GasDiffusionCoeff>"))
            {
                input >> reader;
                IFP_GasDiffCoeff = atof(reader);
                cout << "-----> IFP: Gas Diffusion Coefficient: "
                        << IFP_GasDiffCoeff << "\n";
            }
            if (reader == string("<CO2ResidualSaturation>"))
            {
                input >> reader;
                IFP_CO2ResidSat = atof(reader);
                cout << "-----> IFP: Residual Saturation of CO2: "
                        << IFP_CO2ResidSat << "\n";
            }
            if (reader == string("<MolecularDiffusionCoefficient>"))
            {
                input >> reader;
                IFP_MolecularDiffusionCoefficient = atof(reader);
                cout << "-----> IFP: Molecular diffusion coefficient: "
                        << IFP_MolecularDiffusionCoefficient << "\n";
            }
            if (reader == string("<ViscosityWettingFluid>"))
            {
                input >> reader;
                IFP_ViscosityWettingFluid = atof(reader);
                cout << "-----> IFP: Viscosity of the wetting phase fluid: "
                        << IFP_ViscosityWettingFluid << "\n";
            }
            if (reader == string("<ViscosityNonWettingFluid>"))
            {
                input >> reader;
                IFP_ViscosityNonWettingFluid = atof(reader);
                cout << "-----> IFP: Viscosity of the non wetting phase fluid: "
                        << IFP_ViscosityNonWettingFluid << "\n";
            }

        }
        input.close();
    }

};

class InterfaceProblemProperties
//Class for Interface Problem Properties
//Integrate Parameters for Probabilistic Collocation Method (iPCM).
//Contained design and uncertainty parameters
{
public:
    //Interface Problem Properties (IPP):
    float IPP_DepthBOR; // Depth BOR
    float IPP_InjectionWellRate; // Injection Well Rate
    float IPP_InjectionWindowSize; // Injection Well Window Size
    float IPP_UpperPressure; // Pressure at a top boundary
    float IPP_LowerPressure; // Pressure at a lower boundary
    float IPP_InfiltrationRate; // A Infiltration rate
    float IPP_MaxTimeStepSize; // Maximum time step size
    float IPP_InfiltrationStartTime; // time to stop an infiltration
    float IPP_InfiltrationEndTime; // time to stop an infiltration
    float IPP_DiscretizationLength;

    InterfaceProblemProperties(const char* ipp_filename)
    //Initialization of IPP Parameters
    {
        using namespace std;
        std::cout
                << "-----> IPP: Interface Soil Properties Initialization ...\n";
        //IPP input file defenition
        ifstream input;

        //IPP file check
        input.open(ipp_filename);
        if (!input)
        {
            cout << "\n";
            cout << "-----> IPP: Fatal error! - Data read \n";
            cout << "-----> IPP: Could not open the input data file: \""
                    << ipp_filename << "\n";
        }

        //iPCM input file reading
        char reader[100]; // variable for input value
        //double K;
        while (!input.eof())
        {
            input >> reader;
            //if (reader==string("<BoundaryAndInitialConditions>"))
            //cout << "-----> IPP: Boundary and Initial Conditions reading ... \n";
            //IPP perameters initialization:
            if (reader == string("<DepthBOR>"))
            {
                input >> reader;
                IPP_DepthBOR = atof(reader);
                cout << "-----> IPP: Depth BOR: " << IPP_DepthBOR << "\n";
            }
            if (reader == string("<InjectionWellRate>"))
            {
                input >> reader;
                IPP_InjectionWellRate = atof(reader);
                cout << "-----> IPP: Injection Well Rate: "
                        << IPP_InjectionWellRate << "\n";
            }
            if (reader == string("<InjectionWellWindowSize>"))
            {
                input >> reader;
                IPP_InjectionWindowSize = atof(reader);
                cout << "-----> IPP: Injection Well Window Size: "
                        << IPP_InjectionWindowSize << "\n";
            }
            if (reader == string("<UpperPressure>"))
            {
                input >> reader;
                IPP_UpperPressure = atof(reader);
                cout << "-----> IPP: Upper pressure: "
                        << IPP_UpperPressure << "\n";
            }
            if (reader == string("<LowerPressure>"))
            {
                input >> reader;
                IPP_LowerPressure = atof(reader);
                cout << "-----> IPP: Lower pressure: "
                        << IPP_LowerPressure << "\n";
            }
            if (reader == string("<InfiltrationRate>"))
            {
                input >> reader;
                IPP_InfiltrationRate = atof(reader);
                cout << "-----> IPP: Infiltration rate: "
                        << IPP_InfiltrationRate << "\n";
            }
            if (reader == string("<MaxTimeStepSize>"))
            {
                input >> reader;
                IPP_MaxTimeStepSize = atof(reader);
                cout << "-----> IPP: Maximum time step size: "
                        << IPP_MaxTimeStepSize << "\n";
            }
            if (reader == string("<InfiltrationStartTime>"))
            {
                input >> reader;
                IPP_InfiltrationStartTime = atof(reader);
                cout << "-----> IPP: Start time of infiltration: "
                        << IPP_InfiltrationStartTime << "\n";
            }
            if (reader == string("<InfiltrationEndTime>"))
            {
                input >> reader;
                IPP_InfiltrationEndTime = atof(reader);
                cout << "-----> IPP: End time of infiltration: "
                        << IPP_InfiltrationEndTime << "\n";
            }
            if (reader == string("<DiscretizationLength>"))
            {
                input >> reader;
                IPP_DiscretizationLength = atof(reader);
                cout << "-----> IPP: Discretization length: "
                        << IPP_DiscretizationLength << "\n";
            }
        }
        input.close();
    }

};

} // end namespace
#endif
