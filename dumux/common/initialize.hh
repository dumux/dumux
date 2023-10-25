// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief The initialize function to be called before using Dumux
 */
#ifndef DUMUX_COMMON_INITIALIZE_HH
#define DUMUX_COMMON_INITIALIZE_HH

#include <string>
#include <algorithm>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_TBB
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/global_control.h>

#ifndef DOXYGEN
namespace Dumux::Detail {

class TBBGlobalControl
{
public:
    static oneapi::tbb::global_control& instance(int& argc, char* argv[])
    {
        int maxNumThreads = oneapi::tbb::info::default_concurrency();
        if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
            maxNumThreads = std::max(1, std::stoi(std::string{ dumuxNumThreads }));

        static oneapi::tbb::global_control global_limit(
            oneapi::tbb::global_control::max_allowed_parallelism, maxNumThreads
        );

        return global_limit;
    }
};

} // namespace Dumux::Detail
#endif // DOXYGEN

#endif // HAVE_TBB


#if DUMUX_HAVE_OPENMP
#include <omp.h>
#endif // DUMUX_HAVE_OPENMP


#if DUMUX_HAVE_KOKKOS
#include <Kokkos_Core.hpp>

#ifndef DOXYGEN
namespace Dumux::Detail {

class KokkosScopeGuard
{
public:
    static Kokkos::ScopeGuard& instance(int& argc, char* argv[])
    {
        Kokkos::InitArguments arguments;
        if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
            arguments.num_threads = std::max(1, std::stoi(std::string{ dumuxNumThreads }));

        static Kokkos::ScopeGuard guard(arguments);
        return guard;
    }
};

} // namespace Dumux::Detail
#endif // DOXYGEN

#endif // DUMUX_HAVE_KOKKOS

namespace Dumux {

void initialize(int& argc, char* argv[])
{
    // initialize MPI if available
    // otherwise this will create a sequential (fake) helper
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_TBB
    // initialize TBB and keep global control alive
    Detail::TBBGlobalControl::instance(argc, argv);
#endif

#if DUMUX_HAVE_OPENMP
    if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
        omp_set_num_threads(
            std::max(1, std::stoi(std::string{ dumuxNumThreads }))
        );
#endif

#if DUMUX_HAVE_KOKKOS
    // initialize Kokkos (command line / environmental variable interface)
    Detail::KokkosScopeGuard::instance(argc, argv);
#endif
}

} // end namespace Dumux

#endif
