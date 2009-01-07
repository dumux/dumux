/* ************************************************************ */
/* Dieses Programm berechnet Plots fuer konstitutive            */
/* Beziehungen fuer nichtisotherme 3-Phasen-3-Komponenten-      */
/* Stroemungen                                                  */
/* ************************************************************ */

#include <iostream>
#include <vector>
#include <iomanip>
#include <string>

#include <boost/format.hpp>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>


/* prototypes of fortran routines */
extern "C" double db_(double* Temp, double* pg_MPa);
extern "C" double hb_ (double* Temp, double* rho); 

// C prototype
extern "C" double SolCO2inWater(double Temp, double pg, double X_NaCl);

const char *progname;

void usage()
{
    fprintf(stderr,
            "Usage: %s T1 T2 numT p1 p2 numP [hiresT1 hiresT2 hiresNumT hiresP1 hiresP2 hiresNumP]* > co2tables.h\n"
            "\n"
            "This program generates enthalpy, density, and solubility tables\n"
            "for CO2 depending on temperature and pressure.\n"
            "\n"
            "Parameters are:\n"
            "   T1     Low temperature for the tables [K]\n"
            "   T2     High temperature for the tables [K]\n"
            "   numT   Number of sampling points between T1 and T2\n"
            "   p1     Low pressure for the tables [MPa]\n"
            "   p2     High pressure for the tables [MPa]\n"
            "   numP   Number of sampling points between p1 and p2\n",
            progname);
    exit(1);
};

struct TableEntry
{
    double density;
    double enthalpy;
    double solubility;
};

struct Table : public std::vector<std::vector<TableEntry> >
{
    enum What {
        Density,
        Enthalpy,
        Solubility
    };

    double minP;
    double maxP;
    int    numP;

    double transitionNorth;
    double transitionSouth;
    
    double minT;
    double maxT;
    int    numT;

    double transitionWest;
    double transitionEast;

    int         hiresIndex;
    Table      *hiresTable;

    Table(int argc, const char **argv, int hiresIdx = 0)
        {
            if (argc%6 != 0) {
                usage();
            }
            
            hiresIndex = hiresIdx;
            hiresTable = NULL;
            
            minT = strtod(argv[0], NULL);
            maxT = strtod(argv[1], NULL);
            numT = strtol(argv[2], NULL, 10);
            minP = strtod(argv[3], NULL);
            maxP = strtod(argv[4], NULL);
            numP = strtol(argv[5], NULL, 10);

            transitionNorth = 0.0;
            transitionSouth = 0.0;
            transitionWest = 0.0;
            transitionEast = 0.0;
            
            resize(numP);
            for (int i = 0; i < size(); ++i) {
                operator[](i).resize(numT);
            }

            if (argc > 6) {
                hiresTable = new Table(argc - 6, argv + 6, hiresIndex + 1);
            }
        };

    // make sure the extend of the high-resolution table is adjusted
    // to the boundaries of low-resolution cells.
    void adjustHiresRange()
        {
            if (!hiresTable)
                return;
            
            // round to lower indices for the minumum pressure and
            // temperature.
            int minPIdx = (int) trunc((numP-1)*(hiresTable->minP - minP)/(maxP - minP));
            int minTIdx = (int) trunc((numT-1)*(hiresTable->minT - minT)/(maxT - minT));

            // round to upper indices for the maximum pressure and
            // temperature.
            int maxPIdx = (int) ceil((numP-1)*(hiresTable->maxP - minP)/(maxP - minP));
            int maxTIdx = (int) ceil((numT-1)*(hiresTable->maxT - minT)/(maxT - minT));
            
            // calculate locations of the indices and use them for the
            // highres area
            hiresTable->minP = minP + minPIdx*(maxP - minP)/(numP-1);
            hiresTable->maxP = minP + maxPIdx*(maxP - minP)/(numP-1);

            hiresTable->minT = minT + minTIdx*(maxT - minT)/(numT-1);
            hiresTable->maxT = minT + maxTIdx*(maxT - minT)/(numP-1);

            assert(0 <= minPIdx);
            assert(minPIdx < maxPIdx);
            assert(maxPIdx < numP);

            assert(0 <= minTIdx);
            assert(minTIdx < maxTIdx);
            assert(maxTIdx < numT);
            
            if (minPIdx > 0) {
                hiresTable->transitionWest = 1e6*0.05*(hiresTable->maxP - hiresTable->minP);
            }
            if (maxPIdx < numP - 1) {
                hiresTable->transitionEast = 1e6*0.05*(hiresTable->maxP - hiresTable->minP);
            }
            if (minTIdx > 0) {
                hiresTable->transitionNorth = 0.05*(hiresTable->maxT - hiresTable->minT);
            }
            if (maxTIdx < numT - 1) {
                hiresTable->transitionSouth = 0.05*(hiresTable->maxT - hiresTable->minT);
            }
            
            // if we've got nested hires ranges, do the same
            // recursively
            hiresTable->adjustHiresRange();
        };

    
    // fill the table with the values of the material law
    void calculate()
        {
            std::cerr << "Calculating table @ hires index #" << hiresIndex << "\n";
            for (int i = 0; i < numP; ++i)
            { 
                double pg_MPa = minP + ((double) i)/(numP - 1)*(maxP - minP);
                std::cerr << "  p=" << pg_MPa << " (" << (int) 100*(pg_MPa - minP)/(maxP - minP) << "%)                     \n";
                
                for (int j = 0; j < numT; ++j)
                {	
                    double Temp = minT + ((double) j)/(numT - 1)*(maxT - minT);
                    
                    // call the fortran code from Span and Wagner to calculate
                    // the density and the enthalpy
                    double rho  = db_(&Temp, &pg_MPa);
                    double enth = (2.190963e+01 +  hb_(&Temp, &rho))*1E3 ; 
                    
                    // calculate the solubility of CO2 in water
                    double pg_Pa = pg_MPa * 1.0E6;
                    double solu = SolCO2inWater(Temp, pg_Pa, 0.048);
                    //double solu = SolCO2inWater(Temp, pg_Pa, 0.1);
                    
                    operator[](i)[j].density = rho;
                    operator[](i)[j].enthalpy = enth;
                    operator[](i)[j].solubility = solu;
                    
                    std::cerr << "  T=" << Temp << " (" << (int) 100*(Temp - minT)/(maxT - minT) << "%)                 \r";
                    std::cerr.flush();
                }
            }

            if (hiresTable)
                hiresTable->calculate();
        }

    std::string name(const std::string &fieldName)
        {
            if (hiresIndex == 0) {
                std::string result = "Tabulated";
                result += fieldName;
                return result;
            }

            return (boost::format("HiRes%d%s")%hiresIndex%fieldName).str();
        };

    void printAll()
        {
            if (hiresTable)
                hiresTable->printAll();
            
            print(Density);
            print(Enthalpy);
            print(Solubility);
        };
    
    void print(What what)
        {
            std::string fieldName = "";
            switch (what)
            {
            case Density:
                fieldName = "density"; break;
            case Enthalpy:
                fieldName = "enthalpy"; break;
            case Solubility:
                fieldName = "solubility"; break;
            };

            std::string FieldName = fieldName;
            FieldName[0] = toupper(fieldName[0]);
            
            std::cout.setf(std::ios::scientific);
            std::cout.precision(15);

            std::cout << 
                "struct " << name(FieldName) << "Traits {\n"
                "    typedef double Scalar;\n"
                
                "    static const char  *name;\n"
                "    static const int    numX = " << numP << ";\n"
                "    static const Scalar xMin = " << minP*1e6 << ";\n"
                "    static const Scalar xMax = " << maxP*1e6 << ";\n"
                
                "    static const int    numY = " << numT << ";\n"
                "    static const Scalar yMin = " << minT << ";\n"
                "    static const Scalar yMax = " << maxT << ";\n";

            if (hiresIndex > 0) {
                std::cout << 
                    "    static const Scalar transitionNorth = " << transitionNorth << ";\n"
                    "    static const Scalar transitionSouth = " << transitionSouth << ";\n"
                    "    static const Scalar transitionWest = " << transitionWest << ";\n"
                    "    static const Scalar transitionEast = " << transitionEast << ";\n";
            }
                                
            
            std::string HiResName = "HiResDummy"; 
            if (hiresTable) {
                HiResName = hiresTable->name(FieldName);
            }

            std::cout << 
                "    static const " << HiResName << " hires;\n"
                "    static const Scalar vals[numX][numY];\n"
                "};\n";

            std::cout << "\n";
            std::cout << "const char  *" << name(FieldName) << "Traits::name = \"" << fieldName << "\";\n";
            std::cout << "const " << HiResName << "  " << name(FieldName) << "Traits::hires;\n";
            std::cout << "\n";
            std::cout << "const double " << name(FieldName) << "Traits::vals[" << numP << "][" << numT << "] = \n"
                      << "{\n";
            for (int i = 0; i < numP; ++i) {
                std::cout << "    {";
                for (int j = 0; j < numT; ++j) {
                    if (j % 5 == 0) {
                        std::cout << "\n        ";
                    }
                    
                    double val = 0;
                    switch (what)
                    {
                    case Density:
                        val = (*this)[i][j].density; break;
                    case Enthalpy:
                        val = (*this)[i][j].enthalpy; break;
                    case Solubility:
                        val = (*this)[i][j].solubility; break;
                    };
                    std::cout << std::setw(25) << val;
                    
                    if (j != operator[](i).size() - 1)
                        std::cout << ", ";
                }
                std::cout << "\n    }";
                if (i != size() - 1)
                    std::cout << ",\n";
                else
                    std::cout << "\n";
            }
            std::cout << "};\n"
                "\n";

            if (hiresIndex == 0) {
                std::cout <<
                    "typedef TabulatedMaterial2< " << name(FieldName) << "Traits > " << name(FieldName) << ";\n" <<
                    name(FieldName) << " tabulated" << FieldName << ";\n";
            }
            else {
                std::cout <<
                    "typedef TabulatedMaterial2HiRes< " << name(FieldName) << "Traits > " << name(FieldName) << ";\n";
            }
            std::cout << "\n";
        };
};

int main(int argc, const char **argv)
{
    progname = argv[0];
  
    if (argc < 7) {
        usage();
    }
    
    const double eps = 1e-12;
    
    Table table(argc - 1, argv + 1);
    table.adjustHiresRange();

    std::string cmd = argv[0];
    for (int i = 1; i < argc; ++i) {
        cmd += " ";
        cmd += argv[i];
    }
    
	printf("/* Tables for CO2 fluid properties calculated according to Span and\n"
           " * Wagner (1996).\n"
           " *\n"
           " * THIS AN AUTO-GENERATED FILE! DO NOT EDIT IT!\n"
           " *\n"
           " * Temperature range: %.3lf K to %.3lf K, using %d sampling points\n"
           " * Pressure range: %.3lf MPa to %.3lf MPa, using %d sampling points\n"
           " *\n"
           " * Generated using:\n"
           " *\n"
           " * %s\n"
           " */\n"
           "\n"
           "#ifndef CO2_TABLES_HH\n"
           "#define CO2_TABLES_HH\n"
           "\n"
           "#include <dumux/material/tabulatedmaterial2.hh>\n"
           "#include <dumux/material/tabulatedmaterial2hires.hh>\n"
           "\n"
           "#include <assert.h>\n"
           "\n"
           "namespace Dune {\n"
           "namespace Co2Tables {\n"
           "\n"
           "\n"
           "class HiResDummy {\n"
           "public:\n"
           "\n"
           "    HiResDummy() {};\n"
           "\n"
           "    bool applies(double x, double y) const\n"
           "    { return false; }\n"
           "    \n"
           "    bool hiresWeight(double x, double y) const\n"
           "    { return 0.0; }\n"
           "    \n"
           "    bool at(double x, double y) const\n"
           "    { return 0.0; }\n"
           "};\n"
           "\n"
           ,
           table.minT, table.maxT, table.numT,
           table.minP, table.maxP, table.numP,

           cmd.c_str());
    
    table.calculate();
    
	// print the results to the standard output
    table.printAll();
    
/*
    std::cout << 
        "typedef TabulatedMaterial2<double, densityNumP, densityNumT> DensityTable;\n"
        "static const DensityTable tabulatedDensity(densityMinP, densityMaxP,\n"
        "                                           densityMinT, densityMaxT,\n"
        "                                           densityArray);\n"
        "\n"
        "typedef TabulatedMaterial2<double, enthalpyNumP, enthalpyNumT> EnthalpyTable;\n"
        "static const EnthalpyTable tabulatedEnthalpy(enthalpyMinP, enthalpyMaxP,\n"
        "                                             enthalpyMinT, enthalpyMaxT,\n"
        "                                             enthalpyArray);\n"
        "\n"
        "typedef TabulatedMaterial2<double, solubilityNumP, solubilityNumT> SolubilityTable;\n"
        "static const SolubilityTable tabulatedSolubility(solubilityMinP, solubilityMaxP,\n"
        "                                                 solubilityMinT, solubilityMaxT,\n"
        "                                                 solubilityArray);\n";
*/
    std::cout << "} }\n"
              << "\n"
              << "#endif\n";

    return 0;
}
