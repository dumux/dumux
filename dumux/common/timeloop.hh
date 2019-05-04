// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Common
 * \brief Manages the handling of time dependent problems
 */
#ifndef DUMUX_TIME_LOOP_HH
#define DUMUX_TIME_LOOP_HH

#include <algorithm>
#include <queue>
#include <iomanip>

#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup Common
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
 *
 * \note Time and time step sizes are in units of seconds
 */
template<class Scalar>
class TimeLoopBase
{
public:
    //! Abstract base class needs virtual constructor
    virtual ~TimeLoopBase() {};

    /*!
     * \brief Return the time \f$\mathrm{[s]}\f$ before the time integration.
     * To get the time after the time integration you have to add timeStepSize() to
     * time().
     */
    virtual Scalar time() const = 0;

    /*!
     * \brief Returns the suggested time step length \f$\mathrm{[s]}\f$
     */
    virtual Scalar timeStepSize() const = 0;

    /*!
     * \brief Get the maximum possible time step size \f$\mathrm{[s]}\f$
     */
    virtual Scalar maxTimeStepSize() const = 0;

    /*!
     * \brief Advance to the next time step.
     */
    virtual void advanceTimeStep() = 0;

    /*!
     * \brief Set the current time step size to a given value.
     * \param dt The new value for the time step size \f$\mathrm{[s]}\f$
     */
    virtual void setTimeStepSize(Scalar dt) = 0;

    /*!
     * \brief Returns true if the simulation is finished.
     */
    virtual bool finished() const = 0;
};

/*!
 * \ingroup Common
 * \brief The default time loop for instationary simulations
 */
template <class Scalar>
class TimeLoop : public TimeLoopBase<Scalar>
{
public:
    TimeLoop(Scalar startTime, Scalar dt, Scalar tEnd, bool verbose = true)
    : timer_(false)
    {
        reset(startTime, dt, tEnd, verbose);
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
    }

    /*!
     * \brief Tells the time loop to stop tracking the time.
     * \return the wall clock time (CPU time) spent until now
     */
    double stop()
    {
        return timer_.stop();
    }

    /*!
     * \brief Reset the timer
     */
    void resetTimer()
    {
        timer_.reset();
    }

    /*!
     * \brief Reset the time loop
     */
    void reset(Scalar startTime, Scalar dt, Scalar tEnd, bool verbose = true)
    {
        verbose_ =
            verbose &&
            Dune::MPIHelper::getCollectiveCommunication().rank() == 0;

        time_ = startTime;
        endTime_ = tEnd;

        previousTimeStepSize_ = 0.0;
        userSetMaxTimeStepSize_ = std::numeric_limits<Scalar>::max();
        timeStepIdx_ = 0;
        finished_ = false;
        timeAfterLastTimeStep_ = 0.0;
        timeStepWallClockTime_ = 0.0;

        // ensure that dt is not greater than tEnd-startTime
        setTimeStepSize(dt);

        timer_.stop();
        timer_.reset();
    }

    /*!
     * \brief Advance time step.
     */
    void advanceTimeStep() override
    {
        timeStepIdx_++;
        time_ += timeStepSize_;
        previousTimeStepSize_ = timeStepSize_;

        // compute how long the last time step took
        const auto cpuTime = wallClockTime();
        timeStepWallClockTime_ = cpuTime - timeAfterLastTimeStep_;
        timeAfterLastTimeStep_ = cpuTime;

        // ensure that using current dt we don't exceed tEnd in next time step
        setTimeStepSize(timeStepSize_);
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
    Scalar time() const final
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
    {
        endTime_ = t;
        if (verbose_)
            std::cout << "Set new end time to t = " << t << " seconds." << std::endl;
    }

    /*!
     * \brief Returns the current wall clock time (cpu time) spend in this time loop
     */
    double wallClockTime() const
    { return timer_.elapsed(); }

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
    void setTimeStepSize(Scalar dt) final
    {
        using std::min;
        timeStepSize_ = min(dt, maxTimeStepSize());
    }

    /*!
     * \brief Set the maximum time step size to a given value.
     *
     * \param maxDt The new value for the maximum time step size \f$\mathrm{[s]}\f$
     * \note This also updates the time step size
     */
    void setMaxTimeStepSize(Scalar maxDt)
    {
        userSetMaxTimeStepSize_ = maxDt;
        setTimeStepSize(timeStepSize_);
    }

    /*!
     * \brief Returns the suggested time step length \f$\mathrm{[s]}\f$ so that we
     *        don't miss the beginning of the next episode or cross
     *        the end of the simulation.
     */
    Scalar timeStepSize() const final
    { return timeStepSize_; }

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
    bool finished() const override
    {
        return finished_ || endTime_-time_ < 1e-10*time_;
    }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    {
        return finished() || endTime_-time_-timeStepSize_ < 1e-10*timeStepSize_;
    }

    /*!
     * \brief The current maximum time step size
     * \note This gets aligned on every setTimeStepSize call to end time
     *       and other possible check points
     */
    Scalar maxTimeStepSize() const override
    {
        if (finished())
            return 0.0;

        using std::min; using std::max;
        return min(userSetMaxTimeStepSize_, max<Scalar>(0.0, endTime_ - time_));
    }

    /*!
     * \brief State info on cpu time.
     * \note Always call this after TimeLoop::advanceTimeStep()
     */
    void reportTimeStep() const
    {
        if (verbose_)
        {
            const auto cpuTime = wallClockTime();
            const auto percent = std::round( time_ / endTime_ * 100 );
            std::cout << "[" << std::fixed << std::setw( 3 ) << std::setfill( ' ' )
                      << std::setprecision( 0 )  << percent << "%] "
                      << "Time step " << timeStepIdx_ << " done in "
                      << std::setprecision( 6 ) << timeStepWallClockTime_ << " seconds. "
                      << "Wall clock time: " << std::setprecision( 3 ) << cpuTime
                      << ", time: " << std::setprecision( 5 ) << time_
                      << ", time step size: " << std::setprecision( 8 ) << previousTimeStepSize_
                      << std::endl;
        }
    }

    /*!
     * \brief Print final status and stops tracking the time.
     */
    template< class Communicator = Dune::CollectiveCommunication<typename Dune::MPIHelper::MPICommunicator> >
    void finalize(const Communicator& comm = Dune::MPIHelper::getCollectiveCommunication())
    {
        auto cpuTime = timer_.stop();

        if (verbose_)
        {
            std::cout << "Simulation took " << cpuTime << " seconds on "
                      << comm.size() << " processes.\n";
        }

        if (comm.size() > 1)
            cpuTime = comm.sum(cpuTime);

        if (verbose_)
        {
            std::cout << "The cumulative CPU time was " << cpuTime << " seconds.\n";
        }
    }

    //! If the time loop has verbose output
    bool verbose() const
    { return verbose_; }

    /*
     * @}
     */

protected:
    Dune::Timer timer_;
    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    Scalar previousTimeStepSize_;
    Scalar userSetMaxTimeStepSize_;
    Scalar timeAfterLastTimeStep_, timeStepWallClockTime_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};

/*!
 * \ingroup Common
 * \brief A time loop with a check point mechanism
 */
template <class Scalar>
class CheckPointTimeLoop : public TimeLoop<Scalar>
{
public:
    CheckPointTimeLoop(Scalar startTime, Scalar dt, Scalar tEnd, bool verbose = true)
    : TimeLoop<Scalar>(startTime, dt, tEnd, verbose)
    {
        periodicCheckPoints_ = false;
        deltaPeriodicCheckPoint_ = 0.0;
        lastPeriodicCheckPoint_ = startTime;
        isCheckPoint_ = false;
    }

    /*!
     * \brief Advance time step.
     */
    void advanceTimeStep() override
    {
        const auto dt = this->timeStepSize();
        const auto nextTime = this->time()+dt;

        //! Check point management, TimeLoop::isCheckPoint() has to be called after this!
        // if we reached a periodic check point
        if (periodicCheckPoints_ && Dune::FloatCmp::eq(nextTime, lastPeriodicCheckPoint_ + deltaPeriodicCheckPoint_, 1e-7))
        {
            lastPeriodicCheckPoint_ += deltaPeriodicCheckPoint_;
            isCheckPoint_ = true;
        }

        // or a manually set check point
        else if (!checkPoints_.empty() && Dune::FloatCmp::eq(nextTime, checkPoints_.front(), 1e-7))
        {
            checkPoints_.pop();
            isCheckPoint_ = true;
        }

        // if not reset the check point flag
        else
        {
            isCheckPoint_ = false;
        }

        // advance the time step like in the parent class
        TimeLoop<Scalar>::advanceTimeStep();
    }

    /*!
     * \brief The current maximum time step size
     * \note This gets aligned on every setTimeStepSize call to end time
     *       and other possible check points
     */
    Scalar maxTimeStepSize() const override
    {
        using std::min;
        const auto maxCheckPointDt = computeStepSizeRespectingCheckPoints_();
        const auto maxDtParent = TimeLoop<Scalar>::maxTimeStepSize();
        return min(maxDtParent, maxCheckPointDt);
    }

    /*!
     * \brief Set a periodic check point
     * \note You can query if we are at a time check point with isCheckPoint()
     * \param interval Set a periodic checkout every [interval] seconds
     * \param offset time from which the periodic check points are supposed to start (simulation time)
     *        the first checkpoint will be at time = offset.
     * \note If offset is in the past the first check point will be at the next
     *       periodic check point greater or equal than time
     * \note This also updates the time step size and potentially reduces the time step size to meet the next check point
     */
    void setPeriodicCheckPoint(Scalar interval, Scalar offset = 0.0)
    {
        using std::signbit;
        if (signbit(interval))
            DUNE_THROW(Dune::InvalidStateException, "Interval has to be positive!");

        periodicCheckPoints_ = true;
        deltaPeriodicCheckPoint_ = interval;
        lastPeriodicCheckPoint_ = offset;
        while (lastPeriodicCheckPoint_ + interval - this->time() < 1e-14*interval)
            lastPeriodicCheckPoint_ += interval;

        if (this->verbose())
            std::cout << "Enabled periodic check points every " << interval
                      << " seconds with the next check point at " << lastPeriodicCheckPoint_ + interval << " seconds." << std::endl;

        // check if the current time point is a check point
        if (Dune::FloatCmp::eq(this->time(), lastPeriodicCheckPoint_, 1e-7))
            isCheckPoint_ = true;

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

    //! disable periodic check points
    void disablePeriodicCheckPoints()
    { periodicCheckPoints_ = false; }

    //! remove all check points
    void removeAllCheckPoints()
    {
        periodicCheckPoints_ = false;
        while (!checkPoints_.empty())
            checkPoints_.pop();
    }

    /*!
     * \brief Whether now is a time checkpoint
     * \note has to be called after TimeLoop::advanceTimeStep()
     */
    bool isCheckPoint() const
    { return isCheckPoint_; }

    /*!
     * \brief add a checkpoint to the queue
     * \note checkpoints have to be provided in ascending order
     * \param t the check point (in seconds)
     * \note This also updates the time step size and potentially reduces the time step size to meet the next check point
     */
    void setCheckPoint(Scalar t)
    {
        // set the check point
        setCheckPoint_(t);

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

    /*!
     * \brief add checkpoints to the queue from a vector of time points
     * \note checkpoints have to be provided in ascending order
     * \param checkPoints the vector of check points
     * \note This also updates the time step size and potentially reduces the time step size to meet the next check point
     */
    void setCheckPoint(const std::vector<Scalar>& checkPoints)
    { setCheckPoint(checkPoints.begin(), checkPoints.end()); }

    /*!
     * \brief add checkpoints to the queue from a container from the first iterator to the last iterator
     * \note checkpoints have to be provided in ascending order
     * \param first iterator to the first element to be inserted
     * \param last iterator to the one-after-last element to be inserted
     * \note This also updates the time step size and potentially reduces the time step size to meet the next check point
     */
    template<class ForwardIterator>
    void setCheckPoint(ForwardIterator first, ForwardIterator last)
    {
        // set the check points
        for (; first != last; ++first)
            setCheckPoint_(*first);

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

private:
    //! Adds a check point to the queue
    void setCheckPoint_(Scalar t)
    {
        if (t < this->time())
        {
            if (this->verbose())
                std::cerr << "Couldn't insert checkpoint at t = " << t
                          << " because that's in the past! (current simulation time is " << this->time() << ")" << std::endl;
            return;
        }

        if (!checkPoints_.empty())
        {
            if (t < checkPoints_.back())
            {
                if (this->verbose())
                    std::cerr << "Couldn't insert checkpoint as it is earlier than the last check point in the queue.\n"
                              << "Checkpoints can only be inserted in ascending order." << std::endl;
                return;
            }
        }

        checkPoints_.push(t);
        if (this->verbose())
            std::cout << "Set check point at t = " << t << " seconds." << std::endl;
    }

    /*!
     * \brief Aligns dt to the next check point
     */
    Scalar computeStepSizeRespectingCheckPoints_() const
    {
        using std::min;
        auto maxDt = std::numeric_limits<Scalar>::max();

        if (periodicCheckPoints_)
            maxDt = min(maxDt, lastPeriodicCheckPoint_ + deltaPeriodicCheckPoint_ - this->time());

        if (!checkPoints_.empty())
            maxDt = min(maxDt, checkPoints_.front() - this->time());

        return maxDt;
    }

    bool periodicCheckPoints_;
    Scalar deltaPeriodicCheckPoint_;
    Scalar lastPeriodicCheckPoint_;
    std::queue<Scalar> checkPoints_;
    bool isCheckPoint_;
};

} // end namespace Dumux

#endif
