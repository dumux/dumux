#include <config.h>
#include <iostream>

#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameterparser.hh>
#include <dumux/common/loggingparametertree.hh>

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse the input file into the parameter tree
    Dune::ParameterTree tree;
    Dumux::ParameterParser::parseInputFile(argc, argv, tree, "params.input");

    // attach the tree to the logging tree
    Dumux::LoggingParameterTree params(tree);

    // use some default parameters
    bool DUNE_UNUSED(enableGravity) = params.get<bool>("Problem.EnableGravity", true);

    // use some given parameters
    const auto DUNE_UNUSED(cells) = params.get<std::array<int, 2>>("Grid.Cells", {1, 1});
    const auto DUNE_UNUSED(tEnd) = params.get<double>("TimeLoop.TEnd");

    // check the unused keys
    const auto unused = params.getUnusedKeys();
    if (unused.size() != 1)
        DUNE_THROW(Dune::InvalidStateException, "There should be exactly one unused key!");
    else if (unused[0] != "Grid.Bells")
        DUNE_THROW(Dune::InvalidStateException, "Unused key \"Grid.Bells\" not found!");

    params.reportAll();

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
