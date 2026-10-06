//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parameter
 * \brief Test that Parameters::reset discards the parameters and the record of their use
 */
#include <config.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // first run: set and use A and B
    Parameters::init([](Dune::ParameterTree& params){
        params["Problem.A"] = "1";
        params["Problem.B"] = "2";
    });

    if (getParam<int>("Problem.A") != 1 || getParam<int>("Problem.B") != 2)
        DUNE_THROW(Dune::Exception, "First run does not read the parameters it set");
    if (!getParam<bool>("Problem.EnableGravity"))
        DUNE_THROW(Dune::Exception, "First run does not see the global default of Problem.EnableGravity");

    Parameters::reset();

    if (hasParam("Problem.A") || hasParam("Problem.B"))
        DUNE_THROW(Dune::Exception, "Parameters survive the reset");
    if (getParam<int>("Problem.A", 7) != 7)
        DUNE_THROW(Dune::Exception, "A reset parameter does not fall back to the given default");

    bool globalDefaultsCleared = false;
    try { getParam<bool>("Problem.EnableGravity"); }
    catch (const ParameterException&) { globalDefaultsCleared = true; }
    if (!globalDefaultsCleared)
        DUNE_THROW(Dune::Exception, "Global defaults survive the reset");

    // second run: set B again and C, use only C
    Parameters::init([](Dune::ParameterTree& params){
        params["Problem.B"] = "3";
        params["Problem.C"] = "4";
    });

    if (hasParam("Problem.A"))
        DUNE_THROW(Dune::Exception, "Second run inherits a parameter of the first run");
    if (getParam<int>("Problem.C") != 4)
        DUNE_THROW(Dune::Exception, "Second run does not read the parameters it set");
    if (!getParam<bool>("Problem.EnableGravity"))
        DUNE_THROW(Dune::Exception, "Second run does not see the global defaults");

    // B was used in the first run only, so it has to be reported as unused
    const auto unused = Parameters::getTree().getUnusedKeys();
    if (unused != std::vector<std::string>{"Problem.B"})
    {
        std::ostringstream keys;
        for (const auto& key : unused)
            keys << " " << key;
        DUNE_THROW(Dune::Exception, "Expected Problem.B as the only unused parameter, got:" << keys.str());
    }

    std::ostringstream report;
    Parameters::getTree().reportAll(report);
    if (report.str().find("A = \"1\"") != std::string::npos)
        DUNE_THROW(Dune::Exception, "The report lists a parameter used in the first run:\n" << report.str());

    std::cout << "Parameters::reset test passed" << std::endl;
    return 0;
}
