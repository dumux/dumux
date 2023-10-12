// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Function printing a default usage message
 */
#ifndef DUMUX_DEFAULT_USAGE_MESSAGE_HH
#define DUMUX_DEFAULT_USAGE_MESSAGE_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief Provides a general text block, that is part of error/ help messages.
 *
 * \return The string that is the help / error message.
 */
inline std::string defaultUsageMessage(const std::string& programName)
{
    return  "Usage: " + programName + " [options] \n"
            "Options usually are parameters given to the simulation, \n"
            "and have to be specified with this syntax: \n"
            "\t-GroupName.ParameterName VALUE, for example -TimeLoop.TEnd 100\n"
            "\n"
            "Parameters can also be defined in a parameter file that consists of\n"
            "lines of the form \n"
            "GroupName.ParameterName = VALUE # comment \n"
            "have to be used. More conveniently, group names can be specified in square brackets, \n"
            "such that each following parameter name belongs to that group, \n"
            "[GroupName] \n"
            "ParameterName = VALUE \n"
            "See files named `params.input` in the `test` folder for examples \n"
            "and the Dune documentation of ParameterTreeParser for the format specification. \n"
            "\n"
            "Parameters specified on the command line have priority over those in the parameter file.\n"
            "If no parameter file name is given, './<programname>.input' is chosen as first\n"
            "and './params.input' as second default.\n"
            "\n"
            "Important options include:\n"
            "\t-h, --help                        Print this usage message and exit\n"
            "\t-PrintParameters [true|false]     Print the run-time modifiable parameters _after_ \n"
            "\t                                  the simulation [default: true]\n"
            "\t-ParameterFile FILENAME           File with parameter definitions\n"
            "\t-TimeLoop.Restart RESTARTTIME     Restart simulation from a restart file\n"
            "\n\n";
}

} // end namespace Dumux

#endif
