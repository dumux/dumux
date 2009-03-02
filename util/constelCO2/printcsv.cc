#include <iostream>

#include <dumux/material/constrel/constrelco2.hh>

template <class Traits>
void printArray()
{
    int numT = Traits::numY;
    int numP = Traits::numX;
    double minP = Traits::xMin;
    double maxP = Traits::xMax;
    double minT = Traits::yMin;
    double maxT = Traits::yMax;
    for (int i = 0; i < numT; ++i) {
        for (int j = 0; j < numP; ++j) {
            double T = double(i)/(numT - 1)*(maxT - minT) + minT;
            double p = double(j)/(numP - 1)*(maxP - minP) + minP;
            std::cout << T << " " << p << " " << Traits::vals[j][i] << "\n";
        }
        std::cout << "\n";
    }
};

template <class Traits, class ConstRel>
void printRange(int numP, int numT)
{
    double minP = Traits::xMin;
    double maxP = Traits::xMax;
    double minT = Traits::yMin;
    double maxT = Traits::yMax;

    for (int i = 0; i < numT; ++i) {
        for (int j = 0; j < numP; ++j) {
            double T = double(i)/(numT - 1)*(maxT - minT) + minT;
            double p = double(j)/(numP - 1)*(maxP - minP) + minP;
            std::cout << T << " " << p << " " << ConstRel::enthalpy(T, p) << "\n";
        }
        std::cout << "\n";
    }
}

int main(int argc, char **argv)
{
    typedef Dune::ConstrelCO2 CR;

    if (argc < 2 || argv[1][0] == '0')
        printArray<Dune::Co2Tables::TabulatedEnthalpyTraits>();
    else if (argv[1][0] == '1')
        printArray<Dune::Co2Tables::HiRes1EnthalpyTraits>();
    else if (argv[1][0] == '2')
        printArray<Dune::Co2Tables::HiRes2EnthalpyTraits>();
    /*    else if (argv[1][0] == '3')
          printArray<Dune::Co2Tables::HiRes3EnthalpyTraits>();
    */

    else if (argv[1][0] == 'a')
        printRange<Dune::Co2Tables::TabulatedEnthalpyTraits,
            Dune::ConstrelCO2>(100, 100);

    return 0;
};
