#include <config.h>

#include <iostream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>

namespace Dumux {

Dune::ParameterTree
extractAndTranslateParams(const std::string& paramGroup,
                          const std::vector<std::pair<std::string, std::string>>& translatedKeyList,
                          const std::string& keyPrefix = "",
                          const std::string& translatedKeyPrefix = "")
{
    Dune::ParameterTree translatedParams;

    const std::string emptyString{""};
    for (const auto& [key, translatedKey] : translatedKeyList)
    {
        // if key doesn't exist, empty string is returned
        const auto prefixedKey = keyPrefix == "" ? key : keyPrefix + "." + key;
        auto value = getParamFromGroup<std::string>(paramGroup, prefixedKey, emptyString);

        // if value is not empty string -> translate and insert
        if (value != emptyString)
        {
            const auto prefixedTranslatedKey = translatedKeyPrefix == "" ? translatedKey : translatedKeyPrefix + "." + translatedKey;
            translatedParams[prefixedTranslatedKey] = std::move(value);
        }
    }

    return translatedParams;
}

void checkParameter(const Dune::ParameterTree& params, const std::string& key, const std::string& value)
{
    const auto v = params.get<std::string>(key);
    if (v != value)
        DUNE_THROW(Dune::Exception, "Key: " << key << " returns " << v << " but should return " << value);
}

} // end namespace Dumux

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv, "params_solver.input");

    // the parameter group of our solvers
    std::string paramGroup0{"SubProblem0"};
    std::string paramGroup1{"SubProblem1"};

    // a list of parameters to check and translate
    std::vector<std::pair<std::string, std::string>> solverParamsDumuxToIstl {
        {"Type", "type"},
        {"ResidualReduction", "reduction"},
        {"MaxIterations", "maxIter"},
        {"Preconditioner.Type", "preconditioner.type"}
    };

    // group prefix for linear solvers
    const std::string solverPrefix{"LinearSolver"};

    // create a configuration parameter tree for paramGroup0
    auto configParams0 = extractAndTranslateParams(paramGroup0, solverParamsDumuxToIstl, solverPrefix);
    std::cout << "\ndune-istl solver config for SubProblem0:"
              << "\n-----------------------------------------\n";
    configParams0.report();

    // create a configuration parameter tree for paramGroup1
    auto configParams1 = extractAndTranslateParams(paramGroup1, solverParamsDumuxToIstl, solverPrefix);
    std::cout << "\ndune-istl solver config for SubProblem1:"
              << "\n-----------------------------------------\n";
    configParams1.report();

    // test output
    checkParameter(configParams0, "reduction", "1e-6");
    checkParameter(configParams0, "type", "cgsolver");
    checkParameter(configParams0, "preconditioner.type", "ssor");
    checkParameter(configParams1, "reduction", "1e-8");
    checkParameter(configParams1, "type", "cgsolver");
    checkParameter(configParams1, "preconditioner.type", "amg");
    checkParameter(configParams1, "maxIter", "2000");

    std::cout << "\nConfig correct!" << std::endl;

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e)
{
    std::cout << e << std::endl;
    return 1;
}
