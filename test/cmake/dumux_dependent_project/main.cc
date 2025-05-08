// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/capabilities.hh>

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    const int numProcesses = Dumux::getParam<int>("NumProcesses");

    // do something with dune-common and check MPI
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (const int size = mpiHelper.size(); numProcesses != size)
        DUNE_THROW(Dune::Exception, "Wrong number of processes: have " << size << ", expected " << numProcesses);

    // do something with dune-grid
    static_assert(!Dune::Capabilities::template canCommunicate<void, 0>::v);

    return 0;
}
