//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dumux/io/chrono.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace std::chrono;
    using dseconds = std::chrono::duration<double>;
    using dmilliseconds = std::chrono::duration<double, std::milli>;

    // seconds (test scientific notation and real values)
    if (Chrono::toSeconds("1") != seconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1s") != seconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1.0s") != seconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("42.42s") != dseconds(42.42)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1e3s") != dseconds(1.0e3)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("42.42e3s") != dseconds(42.42e3)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1e-3s") != milliseconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1.5e-3s") != dmilliseconds(1.5)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");

    if (Chrono::toSeconds("1ms") != milliseconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1us") != microseconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1ns") != nanoseconds(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");

    if (Chrono::toSeconds("1min") != minutes(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1h") != hours(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
    if (Chrono::toSeconds("1d") != hours(24)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
#if __cplusplus >= 202002L
    if (Chrono::toSeconds("1y") != years(1)) DUNE_THROW(Dune::InvalidStateException, "Wrong parse result");
#endif

    return 0;
}
