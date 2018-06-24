#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv, "params.input");

    // use some default parameters
    bool enableGravity = getParam<bool>("Problem.EnableGravity", false); // used user default
    if (enableGravity) DUNE_THROW(Dune::InvalidStateException, "Gravity should be false!");

    enableGravity = getParam<bool>("Problem.EnableGravity"); // uses the Dumux default value
    if (!enableGravity) DUNE_THROW(Dune::InvalidStateException, "Gravity should be true!");

    // use some given parameters
    const auto DUNE_UNUSED cells = getParam<std::array<int, 2>>("Grid.Cells", std::array<int, 2>{{1, 1}});
    if (cells[0] != 100 || cells[1] != 100) DUNE_THROW(Dune::InvalidStateException, "Cells should be 100 100!");

    auto tEnd = getParam<double>("TimeLoop.TEnd");
    if (tEnd != 1e6) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e6!");

    tEnd = getParam<double>("TimeLoop.TEnd", 1.0);
    if (tEnd != 1e6) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e6!");

    tEnd = getParamFromGroup<double>("Bulk", "TimeLoop.TEnd", 1.0);
    if (tEnd != 1e5) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e5!");

    tEnd = getParamFromGroup<double>("Bulk", "TimeLoop.TEnd");
    if (tEnd != 1e5) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e5!");

    tEnd = getParamFromGroup<double>("Hulk", "TimeLoop.TEnd", 1.0);
    if (tEnd != 1e6) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e6!");

    tEnd = getParamFromGroup<double>("Hulk", "TimeLoop.TEnd");
    if (tEnd != 1e6) DUNE_THROW(Dune::InvalidStateException, "TEnd should be 1e6!");

    Parameters::print();

    // check the unused keys
    const auto unused = Parameters::getTree().getUnusedKeys();
    if (unused.size() != 1)
        DUNE_THROW(Dune::InvalidStateException, "There should be exactly one unused key!");
    else if (unused[0] != "Grid.Bells")
        DUNE_THROW(Dune::InvalidStateException, "Unused key \"Grid.Bells\" not found!");

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::RangeError &e) {
    std::cout << e << std::endl;
    return 1;
}
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
}
catch (...) {
    std::cout << "Unknown exception!" << std::endl;
    return 1;
}
