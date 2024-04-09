// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Manages the handling of time dependent problems
 */
#ifndef DUMUX_TIME_LOOP_DROPLET_HH
#define DUMUX_TIME_LOOP_DROPLET_HH

#include <algorithm>
#include <queue>
#include <iomanip>
#include <chrono>
#include <type_traits>
#include <initializer_list>
#include <optional>
#include <cmath>

#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>

#include <dune/common/parallel/communication.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/format.hh>
#include <dumux/common/timeloop.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief The default time loop for instationary simulations
 */
template <class Scalar>
class TimeLoopDroplet : public TimeLoop<Scalar>
{
    using ParentType = TimeLoop<Scalar>;
public:
    TimeLoopDroplet(Scalar startTime, Scalar dt, Scalar tEnd, Scalar inkjetPrintingFr, bool verbose = true)
    : ParentType(startTime, dt, tEnd, verbose)
    {
        reset(startTime, dt, tEnd, inkjetPrintingFr);
    }

    template<class Rep1, class Period1,
             class Rep2, class Period2,
             class Rep3, class Period3>
    TimeLoopDroplet(std::chrono::duration<Rep1, Period1> startTime,
             std::chrono::duration<Rep2, Period2> dt,
             std::chrono::duration<Rep3, Period3> tEnd,
             std::chrono::duration<Rep3, Period3> inkjetPrintingFr,
             bool verbose = true)
    : TimeLoopDroplet(
        Detail::TimeLoop::toSeconds<Scalar>(startTime),
        Detail::TimeLoop::toSeconds<Scalar>(dt),
        Detail::TimeLoop::toSeconds<Scalar>(tEnd),
        Detail::TimeLoop::toSeconds<Scalar>(inkjetPrintingFr),
        verbose
    ){}


    /*!
     * \brief Reset the time loop
     */
    template<class Rep1, class Period1,
             class Rep2, class Period2,
             class Rep3, class Period3>
    void reset(std::chrono::duration<Rep1, Period1> startTime,
               std::chrono::duration<Rep2, Period2> dt,
               std::chrono::duration<Rep3, Period3> tEnd,
               std::chrono::duration<Rep3, Period3> inkjetPrintingFr)
    {
        reset(
            Detail::TimeLoop::toSeconds<Scalar>(startTime),
            Detail::TimeLoop::toSeconds<Scalar>(dt),
            Detail::TimeLoop::toSeconds<Scalar>(tEnd),
            Detail::TimeLoop::toSeconds<Scalar>(inkjetPrintingFr));
    }

    /*!
     * \brief Reset the time loop
     */
    void reset(Scalar startTime, Scalar dt, Scalar tEnd, Scalar inkjetPrintingFr)
    {
        inkjetPrintingFrequency_ = inkjetPrintingFr; // Hz
        dropletDispenseTimeInterval_ = 1 / inkjetPrintingFrequency_;
    }

    /*!
     * \brief The current maximum time step size
     * \note This gets aligned on every setTimeStepSize call to end time
     *       and other possible check points
     */
    Scalar maxTimeStepSize() const override
    {
        using std::min; using std::max;

        if (this->finished())
            return 0.0;
        Scalar temp = dropletDispenseTimeInterval_ - (this->time_ - std::round(this->time_ / dropletDispenseTimeInterval_) * dropletDispenseTimeInterval_);
        Scalar dispenseTimeStepSize = temp > dropletDispenseTimeInterval_ + 1e-10 ? (temp - dropletDispenseTimeInterval_):temp;
std::cout<<"--------dispenseTimeStepSize--------"<<dispenseTimeStepSize<<std::endl;
        return min(this->endTime_ - this->time_, dispenseTimeStepSize);
    }

    /*
     * @}
     */

    /*!
     * \brief Return the inkjet printing frequency \f$\mathrm{[1/s]}\f$.
     * inkjetPrintingFrequency().
     */
    Scalar inkjetPrintingFrequency() const
    { return inkjetPrintingFrequency_; }

    Scalar dropletDispenseTimeInterval() const
    { return dropletDispenseTimeInterval_; }

protected:
    static constexpr Scalar baseEps_ = 1e-10;
    Scalar inkjetPrintingFrequency_;
    Scalar dropletDispenseTimeInterval_;
};

} // end namespace Dumux

#endif
