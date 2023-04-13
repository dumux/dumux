// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/python/common/timeloop.hh>

PYBIND11_MODULE(_common, module)
{
    using namespace Dumux;
    using pybind11::operator""_a;

    // maybe initialize MPI and/or multithreading backend
    int argc = 0;
    char **argv = NULL;
    Dumux::initialize(argc, argv);

    // export time loop
    Python::registerTimeLoop<double>(module);

    // initialize the parameters
    Parameters::init();

    // register parameter initializer for command line arguments
    // Usage in Python: initParameters(sys.argv, {"bla": "blubb",})
    module.def("initParameters", [](
        std::vector<std::string>& argv,
        pybind11::dict params,
        const std::string& inputFileName
    ){
        // convert command line args to compatible format
        std::vector<char *> argv_;
        for (auto& arg : argv)
            argv_.push_back(arg.data());

        // convert command line args to parameter tree
        const auto cmdParams = Parameters::parseCommandLine(argv_.size(), argv_.data());

        const auto setParams = [&](auto& paramTree){
            paramTree = cmdParams;
            for (auto [key, value] : params)
                paramTree[pybind11::cast<std::string>(key)] = pybind11::cast<std::string>(pybind11::str(value));
        };

        if (inputFileName.empty())
            Parameters::init(setParams);
        else
            Parameters::init(inputFileName, setParams, /*input overwrites params*/false);

    }, "argv"_a, "params"_a = pybind11::dict{}, "inputFileName"_a = "");

    module.def("getParam", [](const std::string& name){
        return getParam<std::string>(name);
    }, "name"_a);

    module.def("getParam", [](const std::string& name, const std::string& defaultValue){
        return getParam<std::string>(name, defaultValue);
    }, "name"_a, "default"_a);

    module.def("printParameters", &Parameters::print);
}
