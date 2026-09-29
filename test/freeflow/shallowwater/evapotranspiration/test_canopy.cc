// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The canopy partition of a reference evapotranspiration.
 *
 * What is checked here is not the value of the split but that it conserves the demand and that
 * it reaches its two limits: a bare canopy must hand the whole rate to the ground, because that
 * is what makes the winter half of the year the evaporation-only model, and a closed one must
 * keep nearly all of it, because that is what moves the summer demand onto the root zone. A
 * partition that leaked in between would change the annual water balance without changing any
 * parameter that names it.
 */
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/freeflow/shallowwater/evapotranspiration/canopy.hh>

namespace {

void check(bool condition, const std::string& what)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, what);
}

} // end anonymous namespace

int main()
{
    using namespace Dumux::Evapotranspiration;
    using Canopy = DeciduousCanopy<double>;

    // beech in Luxembourg: out from mid-April, bare from early November
    const Canopy beech(5.0, 0.5, 110.0, 305.0, 20.0);

    check(beech.leafAreaIndex(10.0) == 0.0, "a bare canopy has leaf area in midwinter");
    check(beech.leafAreaIndex(350.0) == 0.0, "a bare canopy has leaf area in December");
    check(Dune::FloatCmp::eq(beech.leafAreaIndex(200.0), 5.0),
          "a closed canopy is not at its maximum leaf area in midsummer");

    // the ramps are linear and meet the plateau where they should
    check(beech.leafAreaIndex(110.0) == 0.0, "the canopy has leaf area before leaf-out");
    check(Dune::FloatCmp::eq(beech.leafAreaIndex(120.0), 2.5), "the leaf-out ramp is not linear");
    check(Dune::FloatCmp::eq(beech.leafAreaIndex(130.0), 5.0), "leaf-out does not complete");
    check(Dune::FloatCmp::eq(beech.leafAreaIndex(295.0), 2.5), "the leaf-fall ramp is not linear");
    check(beech.leafAreaIndex(305.0) == 0.0, "the canopy is not bare after leaf fall");

    // the two limits the annual balance turns on
    check(beech.groundFraction(10.0) == 1.0,
          "a bare canopy does not hand the whole demand to the ground");
    check(beech.transpiringFraction(200.0) > 0.9,
          "a closed canopy intercepts less than nine tenths of the demand");

    // The split is Beer's law and not merely something monotonic, pinned against the closed
    // form at three leaf areas: none, half way up the ramp, and closed. A partition that
    // dropped the exponential but kept the right limits would pass every check above.
    check(Dune::FloatCmp::eq(beech.transpiringFraction(200.0), 1.0 - std::exp(-0.5*5.0), 1e-12),
          "a closed canopy does not intercept the Beer's-law share");
    check(Dune::FloatCmp::eq(beech.transpiringFraction(120.0), 1.0 - std::exp(-0.5*2.5), 1e-12),
          "a half-open canopy does not intercept the Beer's-law share");
    check(Dune::FloatCmp::eq(beech.groundFraction(200.0), std::exp(-0.5*5.0), 1e-12),
          "the ground does not receive the transmitted share");

    // and the extinction coefficient has to be the thing in the exponent
    const Canopy thin(5.0, 0.2, 110.0, 305.0, 20.0);
    check(Dune::FloatCmp::eq(thin.transpiringFraction(200.0), 1.0 - std::exp(-0.2*5.0), 1e-12),
          "the extinction coefficient does not enter the exponent");
    check(thin.transpiringFraction(200.0) < beech.transpiringFraction(200.0),
          "a thinner canopy intercepts at least as much as a dense one");

    // With no woody area the two shares are exhaustive, and stay in bounds all year.
    for (int day = 0; day < 366; ++day)
    {
        const auto d = static_cast<double>(day);
        const auto leaves = beech.transpiringFraction(d);
        check(Dune::FloatCmp::eq(leaves + beech.groundFraction(d), 1.0),
              "with no woody area the canopy partition does not sum to one");
        check(leaves >= 0.0 && leaves <= 1.0,
              "the transpiring share leaves the unit interval");
    }

    // Wood keeps intercepting after the leaves fall, so a dormant stand is not bare ground.
    // With woody area, the two shares sum to less than one: what the wood takes is an energy
    // sink that does not become soil water loss, and handing it to either sink would invent ET.
    const Canopy wooded(5.0, 0.5, 110.0, 305.0, 20.0, 0.7);
    check(wooded.transpiringFraction(10.0) == 0.0, "a leafless stand transpires");
    check(Dune::FloatCmp::eq(wooded.groundFraction(10.0), std::exp(-0.5*0.7), 1e-12),
          "the dormant forest floor does not see the woody-area-attenuated rate");
    check(wooded.groundFraction(10.0) < beech.groundFraction(10.0),
          "woody area does not shade the dormant forest floor at all");
    for (int day = 0; day < 366; ++day)
    {
        const auto d = static_cast<double>(day);
        const auto sum = wooded.transpiringFraction(d) + wooded.groundFraction(d);
        check(sum < 1.0 + 1e-12, "a wooded canopy partitions more than the whole demand");
        check(sum > 0.0, "a wooded canopy partitions none of the demand");
    }
    check(Dune::FloatCmp::eq(wooded.transpiringFraction(200.0), beech.transpiringFraction(200.0)),
          "woody area changes what the leaves transpire");

    // no phenology means an evergreen canopy, which never releases the demand
    const Canopy evergreen(5.0, 0.5);
    for (int day = 0; day < 366; ++day)
        check(Dune::FloatCmp::eq(evergreen.leafAreaIndex(static_cast<double>(day)), 5.0),
              "an evergreen canopy sheds its leaves");

    // A season stated backwards must not be taken for an evergreen canopy, which transpires
    // all winter -- most of an annual water balance, from swapping two arguments.
    bool threw = false;
    try { const Canopy backwards(5.0, 0.5, 305.0, 110.0, 20.0); (void)backwards; }
    catch (const Dune::Exception&) { threw = true; }
    check(threw, "a backwards season is accepted as an evergreen canopy");

    threw = false;
    try { const Canopy instant(5.0, 0.5, 110.0, 305.0, 0.0); (void)instant; }
    catch (const Dune::Exception&) { threw = true; }
    check(threw, "a zero-length leaf-out ramp is accepted and divides by zero");

    threw = false;
    try { const Canopy late(5.0, 0.5, 110.0, 400.0, 20.0); (void)late; }
    catch (const Dune::Exception&) { threw = true; }
    check(threw, "leaf fall past the end of the year is accepted and gets cut short");

    // a run that starts mid-year has to see the same season
    check(Dune::FloatCmp::eq(dayOfYear(0.0, 200.0), 199.0),
          "the day of year does not follow the start of the run");
    check(Dune::FloatCmp::eq(beech.leafAreaIndex(dayOfYear(0.0, 200.0)), 5.0),
          "a run started in midsummer begins with a bare canopy");

    std::cout << "canopy partition ok" << std::endl;
    return 0;
}
