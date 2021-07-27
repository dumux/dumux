// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/common/parameters.hh>
#include <dumux/python/common/timeloop.hh>

PYBIND11_MODULE(_common, module)
{
    using namespace Dumux;
    using pybind11::operator""_a;

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
