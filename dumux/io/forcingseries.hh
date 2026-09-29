// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief One column of a measured, timestamped series as a forcing
 */
#ifndef DUMUX_IO_FORCINGSERIES_HH
#define DUMUX_IO_FORCINGSERIES_HH

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A measured forcing, held piecewise constant over its sampling interval
 *
 * The samples are amounts accumulated over an interval (e.g. a rainfall depth) rather than
 * rates, so they are held constant across that interval. Interpolating between them would
 * shift water in time and flatten the intensity peaks.
 *
 * A missing sample (`NaN`) carries the previous value. A missing measurement is not a
 * measurement of zero, and reading it as one would introduce artificial dry spells.
 */
class ForcingSeries
{
public:
    ForcingSeries() = default;

    /*!
     * \param fileName whitespace-separated columns, with a header line naming the columns
     * \param column the name of the column to read
     * \param interval the sampling interval in seconds, over which each sample is held constant
     * \param scale a factor applied to every sample, e.g. 1e-3/interval to convert millimetres per interval to m/s
     */
    ForcingSeries(const std::string& fileName, const std::string& column,
                  double interval, double scale)
    : interval_(interval)
    {
        std::ifstream file(fileName);
        if (!file)
            DUNE_THROW(Dune::IOError, "Cannot open the forcing series " << fileName);
        if (!(interval > 0.0))
            DUNE_THROW(Dune::IOError, "The sampling interval of " << fileName << " must be positive");

        std::string line;
        std::getline(file, line);
        std::istringstream header(line);
        std::vector<std::string> names;
        for (std::string name; header >> name; )
            names.push_back(name);
        const auto found = std::find(names.begin(), names.end(), column);
        if (found == names.end())
            DUNE_THROW(Dune::IOError, "No column " << column << " in " << fileName);
        const auto offset = std::distance(names.begin(), found);

        while (std::getline(file, line))
        {
            std::istringstream values(line);
            std::string value;
            for (int i = 0; i <= offset && values >> value; ++i) {}
            if (value.empty())
                continue;
            values_.push_back(value == "NaN" ? (values_.empty() ? 0.0 : values_.back())
                                             : std::stod(value)*scale);
        }
        std::cout << "Forcing: " << values_.size() << " intervals of " << interval_ << " s from "
                  << fileName << ", column " << column << std::endl;
    }

    //! The forcing at a given time, zero before the start and after the end of the record
    double operator()(const double time) const
    {
        if (values_.empty() || time < 0.0)
            return 0.0;
        const auto index = static_cast<std::size_t>(time/interval_);
        return index < values_.size() ? values_[index] : 0.0;
    }

    /*!
     * \brief The mean rate over a time step
     *
     * Sampling the series at a single time only yields the correct amount while the time step
     * is shorter than the sampling interval. For longer time steps, sampling skips whole
     * intervals, and the skipped amount never enters the model. Averaging over the time step
     * yields exactly the measured amount for any time step size.
     */
    double mean(const double begin, const double end) const
    {
        if (values_.empty() || !(end > begin))
            return (*this)(begin);

        using std::max, std::min;
        const auto from = max(begin, 0.0);
        const auto to = min(end, duration());
        if (!(to > from))
            return 0.0;

        double integral = 0.0;
        auto index = static_cast<std::size_t>(from/interval_);
        for (auto t = from; t < to && index < values_.size(); ++index)
        {
            const auto next = min(to, (index + 1)*interval_);
            integral += values_[index]*(next - t);
            t = next;
        }
        return integral/(end - begin);
    }

    /*!
     * \brief The mean of a function of the forcing over a time step
     *
     * For a rate, `mean` is the right choice, and past the end of the record the rate is zero.
     * An intensive quantity is different in two respects. A response that is nonlinear in the
     * forcing has to be evaluated on the samples and then averaged, not evaluated on the
     * average; by Jensen's inequality, the two differ. And past the end of the record, the
     * quantity is unknown rather than zero, so this averages over the covered part of the time
     * step only, and holds the last value once no part is covered.
     */
    template<class Transform>
    double meanOf(const double begin, const double end, const Transform& transform) const
    {
        if (values_.empty())
            return transform(0.0);

        using std::max, std::min;
        const auto from = max(begin, 0.0);
        const auto to = min(end, duration());
        if (!(to > from))
            return transform(values_[min(static_cast<std::size_t>(max(from, 0.0)/interval_),
                                         values_.size() - 1)]);

        double integral = 0.0;
        auto index = static_cast<std::size_t>(from/interval_);
        for (auto t = from; t < to && index < values_.size(); ++index)
        {
            const auto next = min(to, (index + 1)*interval_);
            integral += transform(values_[index])*(next - t);
            t = next;
        }
        return integral/(to - from);
    }

    bool empty() const { return values_.empty(); }
    std::size_t size() const { return values_.size(); }
    double duration() const { return values_.size()*interval_; }

private:
    std::vector<double> values_;
    double interval_ = 0.0;
};

} // end namespace Dumux

#endif
