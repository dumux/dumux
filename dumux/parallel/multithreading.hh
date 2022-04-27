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
