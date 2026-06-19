// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parallel
 * \brief Helpers to (de)serialize trivially-copyable values into a byte buffer
 *        (e.g. for point-to-point or collective MPI communication)
 */
#ifndef DUMUX_PARALLEL_PACKING_HH
#define DUMUX_PARALLEL_PACKING_HH

#include <vector>
#include <cstring>
#include <type_traits>

namespace Dumux::Detail {

//! Append the raw bytes of a trivially-copyable value to a byte buffer
template<class T>
  requires std::is_trivially_copyable_v<T>
void packValue(std::vector<char>& buf, const T& value)
{
    const auto* p = reinterpret_cast<const char*>(&value);
    buf.insert(buf.end(), p, p + sizeof(T));
}

//! Read a trivially-copyable value from a byte buffer and advance the cursor
template<class T>
  requires std::is_trivially_copyable_v<T>
T unpackValue(const char*& cursor)
{
    T value;
    std::memcpy(&value, cursor, sizeof(T));
    cursor += sizeof(T);
    return value;
}

} // end namespace Dumux::Detail

#endif
