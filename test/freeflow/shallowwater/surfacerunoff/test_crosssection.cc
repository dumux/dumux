// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Properties of the channel cross-sections.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/io/format.hh>

#include <dumux/freeflow/shallowwater/surfacerunoff/crosssection.hh>

namespace {

int failures = 0;

void check(const std::string& what, const double value, const double reference, const double tolerance)
{
    using std::abs, std::max;
    const auto scale = max(1e-30, max(abs(value), abs(reference)));
    if (!(abs(value - reference)/scale < tolerance))
    {
        std::cerr << Dumux::Fmt::format("FAILED {}: {:.10e} != {:.10e}\n", what, value, reference);
        ++failures;
    }
}

void checkTrue(const std::string& what, const bool condition)
{
    if (!condition)
    {
        std::cerr << "FAILED " << what << "\n";
        ++failures;
    }
}

template<class Section>
void checkMonotone(const std::string& name, const Section& s, const double hMax)
{
    bool areaUp = true, conveyanceUp = true;
    double lastA = -1.0, lastK = -1.0;
    for (int i = 0; i <= 400; ++i)
    {
        const auto h = hMax*i/400.0;
        const auto a = s.area(h), k = s.conveyance(h);
        areaUp = areaUp && (a >= lastA - 1e-12);
        conveyanceUp = conveyanceUp && (k >= lastK - 1e-12);
        lastA = a; lastK = k;
    }
    checkTrue(name + ": area is non-decreasing", areaUp);
    checkTrue(name + ": conveyance is non-decreasing", conveyanceUp);
    check(name + ": dry section is empty", s.area(0.0), 0.0, 1e-14);
    check(name + ": dry section conveys nothing", s.conveyance(0.0), 0.0, 1e-14);
    check(name + ": negative depth is empty", s.area(-1.0), 0.0, 1e-14);
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    using namespace Dumux::SurfaceRunoff;
    using std::pow, std::sqrt, std::abs;

    Dumux::initialize(argc, argv);

    // ------------------------------------------------------------------------- rectangle
    {
        const double b = 20.0;
        const RectangularSection<double> rect(b);
        for (const double h : {0.05, 0.5, 2.0, 7.5})
        {
            check(Fmt::format("rectangle area at h = {}", h), rect.area(h), b*h, 1e-14);
            check(Fmt::format("rectangle perimeter at h = {}", h), rect.wettedPerimeter(h), b + 2*h, 1e-14);
            check(Fmt::format("rectangle conveyance at h = {}", h),
                  rect.conveyance(h), b*h*pow(b*h/(b + 2*h), 2.0/3.0), 1e-14);
        }
        checkMonotone("rectangle", rect, 10.0);

        // a wide rectangle approaches sheet flow, which is the law the model uses by default
        const double h = 0.3;
        const RectangularSection<double> wide(1e6);
        check("wide rectangle recovers h^(5/3) per unit width",
              wide.conveyance(h)/1e6, pow(h, 5.0/3.0), 1e-6);
    }

    // ------------------------------------------------------------------------- trapezoid
    {
        const double b = 4.0, t = 12.0, d = 2.0;
        const TrapezoidalSection<double> trap(b, t, d);
        const auto slopeLength = sqrt(1.0 + pow(0.5*(t - b)/d, 2.0));

        check("trapezoid area at half bank", trap.area(1.0), 0.5*(b + 8.0)*1.0, 1e-14);
        check("trapezoid top width at half bank", trap.topWidth(1.0), 8.0, 1e-14);
        check("trapezoid perimeter at half bank", trap.wettedPerimeter(1.0), b + 2.0*slopeLength, 1e-14);
        check("trapezoid area at bankfull", trap.area(d), 0.5*(b + t)*d, 1e-14);
        check("trapezoid bankfull area accessor", trap.bankFullArea(), 0.5*(b + t)*d, 1e-14);

        // above the bank the section stops widening and gains two vertical walls
        check("trapezoid area above bank", trap.area(3.0), 0.5*(b + t)*d + t*1.0, 1e-14);
        check("trapezoid top width above bank", trap.topWidth(3.0), t, 1e-14);
        check("trapezoid perimeter above bank", trap.wettedPerimeter(3.0),
              b + 2.0*d*slopeLength + 2.0*1.0, 1e-14);

        // continuity of every quantity across the bank, which the flux law needs
        const double eps = 1e-7;
        for (const auto& [name, below, above] : std::vector<std::tuple<std::string, double, double>>{
                {"area", trap.area(d - eps), trap.area(d + eps)},
                {"top width", trap.topWidth(d - eps), trap.topWidth(d + eps)},
                {"perimeter", trap.wettedPerimeter(d - eps), trap.wettedPerimeter(d + eps)},
                {"conveyance", trap.conveyance(d - eps), trap.conveyance(d + eps)}})
            check("trapezoid " + name + " is continuous at the bank", below, above, 1e-5);

        checkMonotone("trapezoid", trap, 10.0);

        // the degenerate shapes both occur in real networks
        const TrapezoidalSection<double> rectAsTrap(b, b, d);
        const RectangularSection<double> rect(b);
        check("trapezoid with equal widths is a rectangle",
              rectAsTrap.conveyance(1.3), rect.conveyance(1.3), 1e-14);

        const TrapezoidalSection<double> triangle(0.0, t, d);
        check("triangle area at bankfull", triangle.area(d), 0.5*t*d, 1e-14);
        checkMonotone("triangle", triangle, 10.0);
    }

    // ------------------------------------------------------------------------- tabulated
    {
        const double b = 4.0, t = 12.0, d = 2.0;
        const TrapezoidalSection<double> trap(b, t, d);

        std::vector<double> depth, area, conveyance;
        for (int i = 0; i <= 60; ++i)
        {
            const auto h = 6.0*i/60.0;
            depth.push_back(h);
            area.push_back(trap.area(h));
            conveyance.push_back(trap.conveyance(h));
        }
        const TabulatedSection<double> table(depth, area, conveyance);

        // exact at the knots, and the interpolation error between them is second order
        for (const double h : {0.0, 0.5, 2.0, 4.0, 6.0})
        {
            check(Fmt::format("table reproduces the knot at h = {}", h), table.area(h), trap.area(h), 1e-12);
            check(Fmt::format("table reproduces the conveyance at h = {}", h),
                  table.conveyance(h), trap.conveyance(h), 1e-12);
        }
        // Between the knots the error is the interpolation's own, and near the bed it is not
        // small on a coarse table because the conveyance goes as h^(5/3) there. What has to
        // hold is that refining the table removes it at second order.
        const auto maxError = [&](const int rows)
        {
            std::vector<double> hs, as, ks;
            for (int i = 0; i <= rows; ++i)
            {
                const auto h = 6.0*i/rows;
                hs.push_back(h); as.push_back(trap.area(h)); ks.push_back(trap.conveyance(h));
            }
            const TabulatedSection<double> t(hs, as, ks);
            double worst = 0.0;
            for (int i = 1; i <= 977; ++i)
            {
                const auto h = 6.0*i/977.0;
                worst = std::max(worst, abs(t.conveyance(h) - trap.conveyance(h))/trap.conveyance(6.0));
            }
            return worst;
        };
        // The rate is 5/3 and not 2, and the reason is physical rather than a defect of the
        // interpolation: the conveyance goes as h^(5/3) at the bed, so its second derivative is
        // unbounded there and the first interval of any table carries the largest error.
        const auto coarse = maxError(60), fine = maxError(120), finer = maxError(240);
        check("refining the table converges at the rate the bed allows",
              coarse/fine, pow(2.0, 5.0/3.0), 0.05);
        check("and it keeps that rate", fine/finer, pow(2.0, 5.0/3.0), 0.05);
        checkTrue("a 60-row table is within 1e-3 of the section it was built from", coarse < 1e-3);

        // the perimeter is recovered by inverting the conveyance, which is how a table that
        // carries K rather than P still gives a hydraulic radius
        check("table recovers the wetted perimeter", table.wettedPerimeter(1.0),
              trap.wettedPerimeter(1.0), 5e-3);

        checkMonotone("table", table, 8.0);

        // above the last row the extrapolation stays monotone rather than stopping
        checkTrue("table extrapolates above its last row",
                  table.area(9.0) > table.area(6.0) && table.conveyance(9.0) > table.conveyance(6.0));

        // a conveyance that falls with depth is data the flux law cannot use, so it is refused
        auto falling = conveyance;
        falling[40] = falling[39]*0.5;
        bool threw = false;
        try { TabulatedSection<double> bad(depth, area, falling); }
        catch (const Dune::Exception&) { threw = true; }
        checkTrue("a decreasing conveyance is rejected", threw);
    }

    if (failures > 0)
    {
        std::cerr << "\n" << failures << " check(s) failed\n";
        return 1;
    }

    std::cout << "\nAll checks passed\n";
    return 0;
}
