// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief One column of a regularly sampled series as a forcing
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
 * Each sample is an amount accumulated over the sampling interval (e.g. a rainfall depth) and
 * is held constant across that interval. A missing sample (`NaN`) takes the previous value.
 */
class ForcingSeries
{
public:
    ForcingSeries() = default;

    /*!
     * \param fileName whitespace-separated columns, with a header line naming the columns and
     *        one row per sampling interval, the first starting at time zero
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
     * Integrates the series over the time step, so the amount entering the model is exact for
     * any time step size. The forcing is zero outside the record.
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
     * The transform is evaluated on the samples and then averaged over the part of the time step
     * covered by the record. Outside the record, the nearest sample is held.
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
