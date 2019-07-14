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
/*!
 * \file
 * \brief Test writing and reading sequence container to and from file
 */
#include <config.h>

#include <ostream>
#include <fstream>
#include <algorithm>

#include <dune/common/parallel/mpihelper.hh>
#include <dumux/io/log.hh>

////////////////////////
// the main function
////////////////////////
int main(int argc, char **argv) try {

    using namespace Dumux;

    Dune::MPIHelper::instance(argc, argv);

    Log::progress << "Initializing Dumux..." << Log::end;
    Log::beginTask("Assembling linear system", LogLevel::progress);
    Log::beginTask("Assemble first sub system", LogLevel::progress);
    LogEventMessage printFound("Found intersection!", LogLevel::info, 3);
    for (int i = 0; i < 10; ++i)
        printFound(" (" + std::to_string(i) + ")");
    printFound.summarize();
    Log::progress << "Assembled residual " << 10 << " seconds." << Log::end;
    Log::endTask();
    Log::beginTask("Assemble second sub system", LogLevel::progress);
    Log::progress << "Assembled residual " << 4 << " seconds." << Log::end;
    Log::endTask();
    Log::progress << "Finishes assembly in " << 14 << " seconds." << Log::end;
    Log::endTask();
    Log::warning << "The simulation will be finished soon! Come back to the office!" << Log::end;
    Log::setSleeping();
    Log::info << "This shouldn't be printed!!" << Log::end;
    Log::setSleeping(false);
    std::ofstream logFile("logfile.log");
    Log::setOutputStream(logFile);
    std::string dump("This is a log dumped to file!");
    Log::progress << dump << Log::end;
    Log::setOutputStream(std::cout);
    Log::progress << "Dumux finished." << Log::end;

    if (Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
    {
        std::ifstream readLog("logfile.log");
        std::stringstream buffer;
        buffer << readLog.rdbuf();
        if (buffer.str() != dump + "\n")
            DUNE_THROW(Dune::IOError, "Something went wrong with the file dump!");
    }

    return 0;
}
catch (const Dune::Exception& e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 1;
}
