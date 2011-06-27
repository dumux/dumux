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
#ifndef GW_PARAMETERS_HH
#define GW_PARAMETERS_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>

/** \todo Please doc me! */
namespace Dumux
{

struct Lens
{
    Dune::FieldMatrix<double,2,2> permeability;
    //double permeability;
    Dune::FieldVector<double,2> lowerLeft;
    Dune::FieldVector<double,2> upperRight;
};

struct Source
{
    Dune::FieldVector<double,2> globalPos;
    double q;
    int index;
};

struct BoundaryCondition
{
    std::vector<bool> neumann;
    std::vector<double> value;
    std::vector<double> endPoint;
    int segmentCount;
};

class InterfaceSoilProperties
//Class for Interface Soil Properties
//Integrate Parameters for Probabilistic Collocation Method (iPCM).
//Contained design and uncertainty parameters
{
public:
    //Interface Soil Properties (ISP):
    double permeability;
    double porosity;
    std::vector <Lens> lenses;

    InterfaceSoilProperties(const char* isp_filename)
    //Initialization of ISP Parameters
    {
        using namespace std;
        double viscosity=0.001;
        double density=1000;

        std::cout << std::endl
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
            if (reader == string("<Viscosity>"))
            {
                input >> reader;
                viscosity = atof(reader);
            }
            if (reader == string("<Density>"))
            {
                input >> reader;
                density = atof(reader);
            }



            //cout << "-----> ISP: Soil Properties reading ... \n";
            //ISP parameters initialization:
            if (reader == string("<GlobalPermeability>"))
            {
                input >> reader;
                permeability = atof(reader);
                cout << "-----> Global Permeability: " << permeability << "\n";
            }
            if (reader == string("<Porosity>"))
            {
                input >> reader;
                porosity = atof(reader);
                cout << "-----> Porosity: " << porosity << "\n";
            }

            if (reader == string("<Lens>"))
            {
                Lens templens;
                templens.permeability=0;

                // templens füllen
                while (reader != string("</Lens>"))
                {
                    input >> reader;
                    if (reader == string("<Permeability>"))
                    {
                        input >> reader;
                        templens.permeability[0][0] = atof(reader);
                        templens.permeability[1][1] = templens.permeability[0][0];
                    }
                    if (reader == string("<Left>"))
                    {
                        input >> reader;
                        templens.lowerLeft[0] = atof(reader);
                    }
                    if (reader == string("<Right>"))
                    {
                        input >> reader;
                        templens.upperRight[0] = atof(reader);
                    }
                    if (reader == string("<Bottom>"))
                    {
                        input >> reader;
                        templens.lowerLeft[1] = atof(reader);
                    }
                    if (reader == string("<Top>"))
                    {
                        input >> reader;
                        templens.upperRight[1] = atof(reader);
                    }
                }
                //testen: ist templens gut?
                //an lenses-vector anhängen.
                lenses.push_back(templens);
                cout << "-----> Lens: Left, Bottom, Right, Top: " << templens.lowerLeft
                        << " " << templens.upperRight << ", Permeability; "<<
                        templens.permeability[0][0] << "\n";
            }

        }
        input.close();
        permeability*=(viscosity/(density*9.81));
        for(int i=0; i< lenses.size();i++)
        {
            lenses[i].permeability[0][0]*=(viscosity/(density*9.81));
            lenses[i].permeability[1][1]*=(viscosity/(density*9.81));
        }

    }

};

class InterfaceFluidProperties
//Class for Interface Fluid Properties
//Integrate Parameters for Probabilistic Collocation Method (iPCM).
//Contained design and uncertainty parameters
{
public:
    //Interface Fluid Properties (IFP):
    float viscosity;// Gas Diffusion Coefficient
    float density;// Residual Saturation of CO2

    InterfaceFluidProperties(const char* ifp_filename)
    //Initialization of IFP Parameters
    {
        using namespace std;
        std::cout << std::endl
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
            if (reader == string("<Viscosity>"))
            {
                input >> reader;
                viscosity = atof(reader);
                cout << "-----> Viscosity: "
                        << viscosity << "\n";
            }
            if (reader == string("<Density>"))
            {
                input >> reader;
                density = atof(reader);
                cout << "-----> Density: "
                        << density << "\n";
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
    Dune::FieldVector<int,2> resolution;
    double depth;
    Dune::FieldVector<double,2> size;
    std::vector <Source> sources;
    BoundaryCondition BCondition[4];
    int plotMode;

    InterfaceProblemProperties(const char* ipp_filename)
    //Initialization of IPP Parameters
    {
        using namespace std;
        std::cout
                << "-----> IPP: Interface Problem Properties Initialization ...\n";
        //IPP input file definition
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
            if (reader == string("<XLength>"))
            {
                input >> reader;
                size[0] = atof(reader);
                cout << "-----> X-Length: " << size[0] << "\n";
            }
            if (reader == string("<YLength>"))
            {
                input >> reader;
                size[1] = atof(reader);
                cout << "-----> Y-Length: " << size[1] << "\n";
            }
            if (reader == string("<ZLength>"))
            {
                input >> reader;
                depth = atof(reader);
                cout << "-----> Z-Length: " << depth << "\n";
            }
            if (reader == string("<XResolution>"))
            {
                input >> reader;
                resolution[0] = atoi(reader);
                cout << "-----> X-Resolution: " << resolution[0] << "\n";
            }
            if (reader == string("<YResolution>"))
            {
                input >> reader;
                resolution[1] = atoi(reader);
                cout << "-----> Y-Resolution: " << resolution[1] << "\n";
            }
            if (reader == string("<PlotMode>"))
            {
                input >> reader;
                plotMode = atoi(reader);
                cout << "-----> Plot-mode: " << plotMode << "\n";
            }
            if (reader == string("<Source>"))
            {
                Source tempSource;
                tempSource.q=0;

                while (reader != string("</Source>"))
                {
                    input >> reader;
                    if (reader == string("<X>"))
                    {
                        input >> reader;
                        tempSource.globalPos[0] = atof(reader);
                    }
                    if (reader == string("<Y>"))
                    {
                        input >> reader;
                        tempSource.globalPos[1] = atof(reader);
                    }
                    if (reader == string("<Q>"))
                    {
                        input >> reader;
                        tempSource.q = atof(reader);
                    }
                }
                //an sources-vector anhängen.
                sources.push_back(tempSource);
                cout << "-----> Sink/Source: Position: " << tempSource.globalPos
                        << " " << ", Q = "<< tempSource.q << "\n";
            }

            if (reader == string("<Boundary>"))
            {

                int boundaryIndex = -1;
                while (boundaryIndex == -1)
                {
                    input >> reader;
                    if (reader == string("<Top/>"))
                        boundaryIndex = 0;
                    if (reader == string("<Bottom/>"))
                        boundaryIndex = 1;
                    if (reader == string("<Left/>"))
                        boundaryIndex = 2;
                    if (reader == string("<Right/>"))
                        boundaryIndex = 3;
                }

                while (reader != string("</Boundary>"))
                {
                    input >> reader;

                    if (reader == string("<Top/>"))
                        boundaryIndex = 0;
                    if (reader == string("<Bottom/>"))
                        boundaryIndex = 1;
                    if (reader == string("<Left/>"))
                        boundaryIndex = 2;
                    if (reader == string("<Right/>"))
                        boundaryIndex = 3;

                    if (reader == string("<Type>"))
                    {
                        input >> reader;
                        while (reader != string("</Type>"))
                        {
                            if (reader == string("neumann"))
                                BCondition[boundaryIndex].neumann.push_back(true);
                            else
                                BCondition[boundaryIndex].neumann.push_back(false);
                            input >> reader;
                        }
                    }

                    if (reader == string("<Value>"))
                    {
                        input >> reader;
                        while (reader != string("</Value>"))
                        {
                            BCondition[boundaryIndex].value.push_back(atof(reader));
                            input >> reader;
                        }
                    }

                    if (reader == string("<EndPoint>"))
                    {
                        input >> reader;
                        while (reader != string("</EndPoint>"))
                        {
                            BCondition[boundaryIndex].endPoint.push_back(atof(reader));
                            input >> reader;
                        }
                    }
                }

                for (int i=0; i<4 ; i++)
                {
                    BCondition[i].segmentCount=std::min(std::min(
                            BCondition[i].value.size(),
                            BCondition[i].neumann.size()),
                            BCondition[i].endPoint.size()+1);
                    if (!BCondition[i].segmentCount)
                    {
                        BCondition[i].value.resize(1);
                        BCondition[i].value[0] = 0;
                        BCondition[i].neumann.resize(1);
                        BCondition[i].neumann[0] = true;
                        BCondition[i].endPoint.resize(0);
                        BCondition[i].segmentCount=1;
//                        std::cout << "Invalid Boundary condition in Boundary "<<i<<". Set to no-flow."<<std::endl;
                    }

                    std::cout<< "-----> Boundary Conditions for ";
                    switch(i)
                    {
                    case 0: std::cout<< "top"; break;
                    case 1: std::cout<< "bottom"; break;
                    case 2: std::cout<< "left"; break;
                    case 3: std::cout<< "right";
                    }
                    std::cout << " Boundary:"<<std::endl;
                    for (int j=0; j<BCondition[i].segmentCount ;j++)
                    {
                        if (BCondition[i].neumann[j])
                            std::cout << "        " << "Neumann   ";
                        else
                            std::cout << "        " << "Dirichlet ";
                        std::cout << BCondition[i].value[j] <<" ";
                        if (j < BCondition[i].segmentCount-1)
                            std::cout << BCondition[i].endPoint[j];
                        std::cout << std::endl;
                    }
                }
            }
        }
        input.close();

        // Calculate for each source the containing element
        for (int sourceNumber = 0; sourceNumber != sources.size(); sourceNumber++)
        {
            sources[sourceNumber].index = std::floor(sources[sourceNumber].globalPos[0]
               * resolution[0]/size[0])
               + std::floor(sources[sourceNumber].globalPos[1]*resolution[1]/size[1])
               * resolution[0];
        }
    }

};

} // end namespace
#endif
