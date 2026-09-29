// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Writers for time series, e.g. of a quantity at a probe location, as text files
 */
#ifndef DUMUX_IO_TIMESERIESWRITER_HH
#define DUMUX_IO_TIMESERIESWRITER_HH

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Writes a time series of a single quantity as a two-column text file
 *
 * Every line is flushed, so the file can be read while a simulation is running.
 */
class TimeSeriesWriter
{
public:
    TimeSeriesWriter(const std::string& fileName, const std::string& quantity, const std::string& unit)
    : file_(fileName)
    {
        file_ << "# time[s] " << quantity << "[" << unit << "]\n";
        file_ << std::scientific << std::setprecision(8);
    }

    void write(double time, double value)
    { file_ << time << " " << value << std::endl; }

private:
    std::ofstream file_;
};

/*!
 * \ingroup InputOutput
 * \brief Writes a time series of several quantities as a text file with one column per quantity
 *
 * Writing related quantities on a common time axis, e.g. the terms of a balance, avoids
 * interpolating them onto each other afterwards.
 */
class TableWriter
{
public:
    TableWriter(const std::string& fileName, const std::vector<std::string>& columns)
    : file_(fileName)
    {
        file_ << "# time[s]";
        for (const auto& c : columns)
            file_ << " " << c;
        file_ << "\n" << std::scientific << std::setprecision(8);
    }

    void write(double time, const std::vector<double>& values)
    {
        file_ << time;
        for (const auto v : values)
            file_ << " " << v;
        file_ << std::endl;
    }

private:
    std::ofstream file_;
};

} // end namespace Dumux

#endif
