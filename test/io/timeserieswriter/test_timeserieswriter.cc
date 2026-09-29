// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Test for the time series writers: headers, column order and precision
 */
#include <config.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/io/timeserieswriter.hh>

namespace {

void check(bool condition, const std::string& what)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, "Check failed: " << what);
}

bool close(double a, double b)
{ return std::abs(a - b) <= 1e-8*std::abs(b); }

std::vector<std::string> readLines(const std::string& fileName)
{
    std::ifstream file(fileName);
    std::vector<std::string> lines;
    for (std::string line; std::getline(file, line); )
        lines.push_back(line);
    return lines;
}

std::vector<double> readNumbers(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<double> numbers;
    for (double value; stream >> value; )
        numbers.push_back(value);
    return numbers;
}

} // end anonymous namespace

int main()
{
    using namespace Dumux;

    {
        TimeSeriesWriter writer("test_timeserieswriter_single.dat", "discharge", "m^3/s");
        writer.write(0.0, 0.0);
        writer.write(10.5, 1.23456789e-3);

        // the writer flushes every line, so the file is complete while the writer is alive
        const auto lines = readLines("test_timeserieswriter_single.dat");
        check(lines.size() == 3, "a header and two lines are written");
        check(lines[0] == "# time[s] discharge[m^3/s]", "the header names the quantity and its unit");

        const auto last = readNumbers(lines[2]);
        check(last.size() == 2 && close(last[0], 10.5) && close(last[1], 1.23456789e-3),
              "time and value are written with eight significant digits");
    }

    {
        TableWriter writer("test_timeserieswriter_table.dat", {"rain[m^3]", "outflow[m^3]", "storage[m^3]"});
        writer.write(60.0, {1.0, 0.25, 0.75});

        const auto lines = readLines("test_timeserieswriter_table.dat");
        check(lines.size() == 2, "a header and one line are written");
        check(lines[0] == "# time[s] rain[m^3] outflow[m^3] storage[m^3]", "the header lists the columns in order");

        const auto values = readNumbers(lines[1]);
        check(values.size() == 4 && close(values[0], 60.0) && close(values[1], 1.0)
              && close(values[2], 0.25) && close(values[3], 0.75),
              "the values follow the time in column order");
    }

    std::cout << "All time series writer checks passed" << std::endl;
    return 0;
}
