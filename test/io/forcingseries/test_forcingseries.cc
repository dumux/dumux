// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Averaging a measured forcing over a time step.
 *
 * The two averages in this class answer different questions and must not be swapped. A rate is
 * averaged so that a long step puts exactly the measured depth in, and running past the end of
 * the record contributes nothing because there is no more rain to add. An intensive reading
 * driving a non-linear response is averaged *after* the response is evaluated, and past the end
 * of the record it is unknown rather than zero. Getting either backwards changes a water balance
 * without changing any parameter that names it.
 */
#include <config.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/io/forcingseries.hh>

namespace {

void check(bool condition, const std::string& what)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, what);
}

std::string writeSeries(const std::string& name, const std::string& body)
{
    std::ofstream file(name);
    file << "datetime value\n" << body;
    return name;
}

} // end anonymous namespace

int main()
{
    using Dumux::ForcingSeries;

    // four intervals of 100 s
    const auto path = writeSeries("test_forcingseries.dat",
                                  "t0 1.0\nt1 3.0\nt2 2.0\nt3 4.0\n");
    const ForcingSeries series(path, "value", 100.0, 1.0);
    check(series.size() == 4, "the series did not read four intervals");
    check(Dune::FloatCmp::eq(series.duration(), 400.0), "the series has the wrong duration");

    // a rate: averaging over a step must conserve the total
    check(Dune::FloatCmp::eq(series.mean(0.0, 400.0), 2.5), "the mean over the record is wrong");
    check(Dune::FloatCmp::eq(series.mean(0.0, 200.0), 2.0), "the mean over two intervals is wrong");
    check(Dune::FloatCmp::eq(series.mean(50.0, 150.0), 2.0), "a straddling step averages wrongly");

    // ...and past the end it dilutes, because there is no more of it to add
    check(Dune::FloatCmp::eq(series.mean(0.0, 800.0), 1.25),
          "a rate past the end of the record does not contribute zero");

    // An intensive reading is different in both respects. Jensen: for a convex response the
    // average of the response exceeds the response of the average, always in the same direction.
    const auto closure = [](const double d) { return 1.0/(1.0 + d); };
    const auto averaged = series.meanOf(0.0, 400.0, closure);
    const auto ofAverage = closure(series.mean(0.0, 400.0));
    check(averaged > ofAverage,
          "averaging a convex response does not exceed the response of the average");
    const auto expected = (closure(1.0) + closure(3.0) + closure(2.0) + closure(4.0))/4.0;
    check(Dune::FloatCmp::eq(averaged, expected, 1e-12),
          "the averaged response is not the mean of the per-interval responses");

    // a step inside one interval is just that interval's response
    check(Dune::FloatCmp::eq(series.meanOf(10.0, 90.0, closure), closure(1.0), 1e-12),
          "a step inside one interval does not take that interval's value");

    // Past the end it must hold the last reading, not decay toward zero: for a stomatal closure
    // a deficit of zero means saturated air and fully open stomata, so diluting would silently
    // remove the limit precisely where the record runs out.
    check(Dune::FloatCmp::eq(series.meanOf(400.0, 500.0, closure), closure(4.0), 1e-12),
          "past the end of the record the reading is not held");
    check(Dune::FloatCmp::eq(series.meanOf(1000.0, 1100.0, closure), closure(4.0), 1e-12),
          "far past the end of the record the reading is not held");
    const auto straddling = series.meanOf(350.0, 450.0, closure);
    check(Dune::FloatCmp::eq(straddling, closure(4.0), 1e-12),
          "a step straddling the end of the record is diluted");

    // a gap carries the previous reading rather than reading as a measurement of zero
    const auto gappy = writeSeries("test_forcingseries_gap.dat",
                                   "t0 1.0\nt1 NaN\nt2 2.0\n");
    const ForcingSeries withGap(gappy, "value", 100.0, 1.0);
    check(Dune::FloatCmp::eq(withGap.mean(100.0, 200.0), 1.0),
          "a gap in the record reads as zero instead of carrying the previous value");

    // an empty series must not be indexed
    const ForcingSeries none;
    check(none.empty(), "a default-constructed series is not empty");
    check(Dune::FloatCmp::eq(none.meanOf(0.0, 100.0, closure), closure(0.0)),
          "an empty series does not fall back to a zero reading");

    std::remove(path.c_str());
    std::remove(gappy.c_str());
    std::cout << "forcing series ok" << std::endl;
    return 0;
}
