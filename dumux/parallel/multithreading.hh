// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup Parallel
 * \brief Multithreading in Dumux
 */
#ifndef DUMUX_PARALLEL_MULTITHREADING_HH
#define DUMUX_PARALLEL_MULTITHREADING_HH

#include <type_traits>

#ifndef DUMUX_MULTITHREADING_BACKEND
#define DUMUX_MULTITHREADING_BACKEND Serial
#endif

namespace Dumux::Detail::Multithreading {

namespace ExecutionBackends {

struct Serial {};
struct Cpp {};
struct TBB {};
struct Kokkos {};
struct OpenMP {};

} // end namespace ExecutionBackends

// set the execution backend type
using ExecutionBackend = ExecutionBackends::DUMUX_MULTITHREADING_BACKEND;

} // end namespace Dumux::Detail::Multithreading

namespace Dumux::Multithreading {

/*!
 * \ingroup Parallel
 * \brief Checking whether the backend is serial
 */
inline constexpr bool isSerial()
{
    using namespace Dumux::Detail::Multithreading;
    return std::is_same_v<ExecutionBackends::Serial, ExecutionBackend>;
};

} // end namespace Dumux

#endif
