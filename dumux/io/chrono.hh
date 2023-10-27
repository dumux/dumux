// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Helper functions for working with std::chrono.
 */
#ifndef DUMUX_IO_CHRONO_HH
#define DUMUX_IO_CHRONO_HH

#include <array>
#include <cctype>
#include <chrono>
#include <string>
#include <string_view>
#include <algorithm>

#include <dune/common/exceptions.hh>

namespace Dumux::Chrono {

//! Try to construct a std::chrono::duration from a string
template<typename Rep, typename Period>
void toDuration(std::chrono::duration<Rep, Period>& duration, const std::string& s)
{
    using std::chrono::duration_cast;
    using namespace std::chrono_literals;
    using S = std::chrono::duration<Rep>;
    constexpr std::array<std::pair<std::string_view, S>, 9> unitMap{{
        {"", S(1)}, // assume seconds if no unit is given
        {"s", S(1)},
        {"ns", duration_cast<S>(std::chrono::nanoseconds(1))},
        {"us", duration_cast<S>(std::chrono::microseconds(1))},
        {"ms", duration_cast<S>(std::chrono::milliseconds(1))},
        {"min", duration_cast<S>(std::chrono::minutes(1))},
        {"h", duration_cast<S>(std::chrono::hours(1))},
        // After requiring cpp20, we can use the aliases in std::chrono
        {"d", duration_cast<S>(std::chrono::hours(24))},
        {"y", duration_cast<S>(std::chrono::seconds(31556952))}  // to match cpp20, see https://en.cppreference.com/w/cpp/chrono/duration
    }};

    const auto unitIt = std::find_if(s.rbegin(), s.rend(), [] (const auto& c) { return !std::isalpha(c); });
    const auto unitPos = s.size() - std::distance(s.rbegin(), unitIt);
    const auto number = std::stod(s.substr(0, unitPos));
    const auto unit = s.substr(unitPos);
    const auto conversion = [&] () {
        const auto it = std::find_if(unitMap.begin(), unitMap.end(), [u=std::string_view{unit}] (const auto& p) {
            return p.first == u;
        });
        if (it == unitMap.end())
            DUNE_THROW(Dune::InvalidStateException, "Unsupported unit " << unit << ".");
        return it->second;
    } ();
    duration = duration_cast<std::chrono::duration<Rep, Period>>(number*conversion);
}

//! Try to construct an instance of std::chrono::seconds from a string including a unit suffix
template<typename Rep = double>
std::chrono::duration<Rep> toSeconds(const std::string& s)
{
    std::chrono::duration<Rep> result;
    toDuration(result, s);
    return result;
}

} // end namespace Dumux::Chrono

#endif
