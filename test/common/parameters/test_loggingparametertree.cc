#include <config.h>
#include <iostream>

#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/parameterparser.hh>

namespace Dumux {

namespace Properties
{
NEW_TYPE_TAG(Bla);
NEW_PROP_TAG(TimeLoopTLEnd);
NEW_PROP_TAG(Scalar);
SET_TYPE_PROP(Bla, Scalar, double);
SET_SCALAR_PROP(Bla, TimeLoopTLEnd, 2.0);
}

}

int main (int argc, char *argv[]) try
{
    using namespace Dumux;
    using TypeTag = TTAG(Bla);

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse the input file into the parameter tree
    // Dune::ParameterTree tree;
    // Dumux::ParameterParser::parseInputFile(argc, argv, tree, "params.input");
    //
    // // attach the tree to the logging tree
    // Dumux::LoggingParameterTree params(tree);
    //
    // // use some default parameters
    // bool DUNE_UNUSED(enableGravity) = params.get<bool>("Problem.EnableGravity", true);
    //
    // // use some given parameters
    // const auto DUNE_UNUSED(cells) = params.get<std::array<int, 2>>("Grid.Cells", {1, 1});
    // const auto DUNE_UNUSED(tEnd) = params.get<double>("TimeLoop.TEnd");
    //
    // // check the unused keys
    // const auto unused = params.getUnusedKeys();
    // if (unused.size() != 1)
    //     DUNE_THROW(Dune::InvalidStateException, "There should be exactly one unused key!");
    // else if (unused[0] != "Grid.Bells")
    //     DUNE_THROW(Dune::InvalidStateException, "Unused key \"Grid.Bells\" not found!");
    //
    // params.reportAll();

    Parameters::init(argc, argv, "params.input");

    // use some default parameters
    bool enableGravity = getParam<bool>("Problem.EnableGravity", false); // used user default
    std::cout << enableGravity << std::endl;
    enableGravity = getParam<bool>("Problem.EnableGravity"); // uses the Dumux default value
    std::cout << enableGravity << std::endl;

    // use some given parameters
    const auto DUNE_UNUSED(cells) = getParam<std::array<int, 2>>("Grid.Cells", std::array<int, 2>{{1, 1}});
    auto tEnd = getParam<double>("TimeLoop.TEnd");
    tEnd = getParam<double>("TimeLoop.TEnd", 1.0);
    // const auto tEnd = getParamFromGroup<double>("Bulk", "TimeLoop.TEnd", 1.0, Parameters::simpleLookup);
    tEnd = getParamFromGroup<double>("Bulk", "TimeLoop.TEnd", 1.0, ParamLookup::tree);
    // tEnd = GET_PARAM(TypeTag, double, TimeLoop.TEnd);
    tEnd = GET_PARAM_FROM_GROUP(TypeTag, double, TimeLoop, TLEnd);
    std::cout << tEnd << std::endl;
    tEnd = GET_RUNTIME_PARAM(TypeTag, double, TimeLoop.TEnd);
    tEnd = GET_RUNTIME_PARAM_CSTRING(TypeTag, double, "TimeLoop.TEnd");
    tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeLoop, TEnd);
    tEnd = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, double, "TimeLoop", TEnd);
    std::cout << tEnd << std::endl;

    reportParams();

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
