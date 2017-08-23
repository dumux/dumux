// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Manages the handling of time dependent problems
 */
#ifndef DUMUX_TIME_LOOP_HH
#define DUMUX_TIME_LOOP_HH

#include <algorithm>

#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>

#include "propertysystem.hh"
#include "parameters.hh"

namespace Dumux
{
/*!
 * \ingroup SimControl
 * \brief Manages the handling of time dependent problems.
 *
 * This class facilitates the time management of the simulation.
 * It doesn't manage any user data, but keeps track of what the
 * current time, time step size and "episode" of the
 * simulation is. It triggers the initialization of the problem and
 * is responsible for the time control of a simulation run.
 *
 * The time manager allows to specify a sequence of "episodes" which
 * determine the boundary conditions of a problem. This approach
 * is handy if the problem is not static, i.e. the boundary
 * conditions change over time.
 *
 * An episode is a span of simulated time in which
 * the problem behaves in a specific way. It is characterized by
 * the (simulation) time it starts, its length and a consecutive
 * index starting at 0.
 */
template <class Scalar>
class TimeLoop
{
public:
    TimeLoop(Scalar startTime, Scalar dt, Scalar tEnd, bool verbose = true) : timer_(false)
    {
        verbose_ =
            verbose &&
            Dune::MPIHelper::getCollectiveCommunication().rank() == 0;

        time_ = startTime;
        endTime_ = tEnd;

        timeStepSize_ = dt;
        previousTimeStepSize_ = timeStepSize_;
        maxTimeStepSize_ = std::numeric_limits<Scalar>::max();
        timeStepIdx_ = 0;
        finished_ = false;
    }

    /*!
     *  \name Simulated time and time step management
     * @{
     */

    /*!
     * \brief Tells the time loop to start tracking the time.
     */
    void start()
    {
        timer_.start();
        cpuTime_ = 0.0;
    }

    /*!
     * \brief State info on cpu time.
     */
    void reportTimeStep()
    {
        auto timeStepCpuTime = timer_.elapsed();
        cpuTime_ += timeStepCpuTime;

        if (verbose_)
        {
                std::cout << "Time step " << timeStepIndex() << " done in "
                          << timeStepCpuTime << " seconds. "
                          << "Wall time: " << cpuTime_
                          << ", time: " << time()
                          << ", time step size: " << timeStepSize()
                          << std::endl;
        }

        timer_.reset();
    }

    /*!
     * \brief Advance time step.
     */
    void advanceTimeStep()
    {
        timeStepIdx_++;
        time_ += timeStepSize_;
    }

    /*!
     * \brief Print final status and stops tracking the time.
     */
    void finalize()
    {
        timer_.stop();
        cpuTime_ += timer_.elapsed();

        if (verbose_)
        {
            const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
            std::cout << "Simulation took " << timer_.elapsed() << " seconds on "
                      << comm.size() << " processes.\n"
                      << "The cumulative CPU time was " << comm.sum(cpuTime_) << " seconds.\n";
        }
    }

    /*!
     * \brief Set the current simulated time, don't change the current
     *        time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    void setTime(Scalar t)
    { time_ = t; }

    /*!
     * \brief Set the current simulated time and the time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     * \param stepIdx The new time step index
     */
    void setTime(Scalar t, int stepIdx)
    { time_ = t; timeStepIdx_ = stepIdx; }

    /*!
     * \brief Return the time \f$\mathrm{[s]}\f$ before the time integration.
     * To get the time after the time integration you have to add timeStepSize() to
     * time().
     */
    Scalar time() const
    { return time_; }

    /*!
     * \brief Returns the number of (simulated) seconds which the simulation runs.
     */
    Scalar endTime() const
    { return endTime_; }

    /*!
     * \brief Set the time of simulated seconds at which the simulation runs.
     *
     * \param t The time \f$\mathrm{[s]}\f$ at which the simulation is finished
     */
    void setEndTime(Scalar t)
    { endTime_ = t; }

    /*!
     * \brief Returns the current wall time (cpu time).
     */
    double wallTime() const
    {  return cpuTime_; }

    /*!
     * \brief Set the current time step size to a given value.
     *
     * If the step size would exceed the length of the current
     * episode, the timeStep() method will take care that the step
     * size won't exceed the episode or the end of the simulation,
     * though.
     *
     * \param dt The new value for the time step size \f$\mathrm{[s]}\f$
     */
    void setTimeStepSize(Scalar dt)
    {
        using std::min;
        timeStepSize_ = min(dt, maxTimeStepSize());
    }

    /*!
     * \brief Set the maximum time step size to a given value.
     *
     * \param dt The new value for the maximum time step size \f$\mathrm{[s]}\f$
     */
    void setMaxTimeStepSize(Scalar maxDt)
    { maxTimeStepSize_ = maxDt; }

    /*!
     * \brief Returns the suggested time step length \f$\mathrm{[s]}\f$ so that we
     *        don't miss the beginning of the next episode or cross
     *        the end of the simulation.
     */
    Scalar timeStepSize() const
    { return timeStepSize_; }

    /*!
     * \brief Returns the size of the previous time step \f$\mathrm{[s]}\f$.
     */
    Scalar previousTimeStepSize() const
    { return previousTimeStepSize_; }

    /*!
     * \brief Returns number of time steps which have been
     *        executed since the beginning of the simulation.
     */
    int timeStepIndex() const
    { return timeStepIdx_; }

    /*!
     * \brief Specify whether the simulation is finished
     *
     * \param finished If true the simulation is considered finished
     *                 before the end time is reached, else it is only
     *                 considered finished if the end time is reached.
     */
    void setFinished(bool finished = true)
    { finished_ = finished; }

    /*!
     * \brief Returns true if the simulation is finished.
     *
     * This is the case if either setFinished(true) has been called or
     * if the end time is reached.
     */
    bool finished() const
    { return finished_ || time() >= endTime(); }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    { return finished_ || time() + timeStepSize() >= endTime(); }

    /*!
     * \brief Aligns dt to the episode boundary or the end time of the
     *        simulation.
     */
    Scalar maxTimeStepSize() const
    {
        if (finished())
            return 0.0;

        using std::max;
        using std::min;
        // return min(min(episodeMaxTimeStepSize(),
        //                problem_->maxTimeStepSize()),
        //            max<Scalar>(0.0, endTime() - time()));
        return min(maxTimeStepSize_,
                   max<Scalar>(0.0, endTime() - time()));
    }

    /*
     * @}
     */

private:
    Dune::Timer timer_;
    Scalar time_;
    Scalar endTime_;
    double cpuTime_;

    Scalar timeStepSize_;
    Scalar previousTimeStepSize_;
    Scalar maxTimeStepSize_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};
}

#endif
