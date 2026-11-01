//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parameter
 * \brief Test the log of parameters that fall back to the default value given at the call site
 */
#include <config.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dumux/common/parameters.hh>

namespace {

const std::string fallbackHeader = "# Fallback default parameters used";

// the report section starting with the given header, empty if the report has no such section
std::string reportSection(const std::string& header)
{
    std::ostringstream stream;
    Dumux::Parameters::getTree().reportAll(stream);
    const auto report = stream.str();
    const auto begin = report.find(header);
    if (begin == std::string::npos)
        return {};
    return report.substr(begin, report.find("\n#", begin) - begin);
}

void checkContains(const std::string& section, const std::string& text, bool expected = true)
{
    if ((section.find(text) != std::string::npos) != expected)
        DUNE_THROW(Dune::Exception, "Expected the report section " << (expected ? "" : "not ")
                                    << "to contain '" << text << "':\n" << section);
}

} // end anonymous namespace

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // without Parameters.LogFallbackDefaults the fallback defaults are not logged
    Parameters::init([](Dune::ParameterTree&){});
    getParam<double>("Newton.MaxRelativeShift", 1e-8);
    if (!reportSection(fallbackHeader).empty())
        DUNE_THROW(Dune::Exception, "Fallback defaults are reported although their log is disabled");

    // enable the log in the input file
    Parameters::reset();
    Parameters::init(argc, argv, "params_fallbackdefaults.input");

    getParam<double>("Newton.MaxRelativeShift", 1e-8);
    getParam<double>("Fallback.Sum", 0.1 + 0.2);
    getParam<bool>("Fallback.Bool", true);
    getParam<int>("Fallback.Int", -3);
    getParam<std::string>("Fallback.String", "a b");
    getParam("Fallback.DefaultTypeString", "c");
    getParam<std::array<int, 2>>("Fallback.Array", std::array<int, 2>{{1, 2}});
    getParam<std::vector<double>>("Fallback.Vector", std::vector<double>{0.5, 1e-3});
    getParam<Dune::FieldVector<double, 2>>("Fallback.FieldVector", Dune::FieldVector<double, 2>{0.1, 2.0});
    getParamFromGroup<int>("Group", "Fallback.InGroup", 5);
    Parameters::getTree().get("Fallback.CString", "d");
    getParam<std::filesystem::path>("Fallback.Path", std::filesystem::path("dir/file"));
    getParam<std::array<int, 2>>("Grid.Cells", std::array<int, 2>{{1, 1}});

    // the logged values are read back as the values used
    const auto fallbackSection = reportSection(fallbackHeader);
    std::istringstream fallbackStream(fallbackSection);
    Dune::ParameterTree fallbacks;
    Dune::ParameterTreeParser::readINITree(fallbackStream, fallbacks);

    if (fallbacks.get<double>("Newton.MaxRelativeShift") != 1e-8
        || fallbacks.get<double>("Fallback.Sum") != 0.1 + 0.2
        || fallbacks.get<bool>("Fallback.Bool") != true
        || fallbacks.get<int>("Fallback.Int") != -3
        || fallbacks["Fallback.String"] != "a b"
        || fallbacks["Fallback.DefaultTypeString"] != "c"
        || fallbacks.get<std::array<int, 2>>("Fallback.Array") != std::array<int, 2>{{1, 2}}
        || fallbacks.get<std::vector<double>>("Fallback.Vector") != std::vector<double>{0.5, 1e-3}
        || fallbacks.get<Dune::FieldVector<double, 2>>("Fallback.FieldVector") != Dune::FieldVector<double, 2>{0.1, 2.0}
        || fallbacks.get<int>("Fallback.InGroup") != 5
        || fallbacks["Fallback.CString"] != "d")
        DUNE_THROW(Dune::Exception, "Fallback defaults are not logged with the values used:\n" << fallbackSection);

    checkContains(fallbackSection, "Sum = \"0.30000000000000004\"");
    checkContains(fallbackSection, "MaxRelativeShift = \"1e-08\"");
    checkContains(fallbackSection, "dir/file");

    // parameters set at runtime are not fallback defaults
    if (fallbacks.hasKey("Grid.Cells") || fallbacks.hasKey("Parameters.LogFallbackDefaults"))
        DUNE_THROW(Dune::Exception, "Runtime parameters are logged as fallback defaults:\n" << fallbackSection);
    checkContains(reportSection("# Runtime-specified parameters used"), "LogFallbackDefaults = \"true\"");

    // a reset discards the log, another init can enable it again
    Parameters::reset();
    Parameters::init([](Dune::ParameterTree& params){
        params["Parameters.LogFallbackDefaults"] = "true";
    });
    getParam<int>("Newton.MaxSteps", 20);
    checkContains(reportSection(fallbackHeader), "MaxSteps = \"20\"");
    checkContains(reportSection(fallbackHeader), "[ Fallback ]", false);

    // and disable it again
    Parameters::reset();
    Parameters::init([](Dune::ParameterTree&){});
    getParam<int>("Newton.MaxSteps", 18);
    if (!reportSection(fallbackHeader).empty())
        DUNE_THROW(Dune::Exception, "Fallback defaults are reported although their log is disabled");

    std::cout << "Fallback default parameter test passed" << std::endl;
    return 0;
}
