#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    auto defaultParams = [] (Dune::ParameterTree& p) {
        p["Grid.UpperRight"] = "2 2";
        p["Grid.Cells"] = "20 20";
    };

    // initialize parameter tree
    Parameters::init(argc, argv, defaultParams);
std::cout << getParam<Dune::FieldVector<double, 2>>("Grid.UpperRight") << std::endl;
    Parameters::print();

    return 0;
}
