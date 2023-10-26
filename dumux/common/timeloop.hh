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
#ifndef DUMUX_TIME_LOOP_HH
#define DUMUX_TIME_LOOP_HH

#include <algorithm>
#include <queue>
#include <iomanip>
#include <chrono>
#include <type_traits>
#include <initializer_list>
#include <optional>

#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>

#include <dune/common/parallel/communication.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/format.hh>

namespace Dumux {


#ifndef DOXYGEN
namespace Detail::TimeLoop {

template<typename Scalar, typename R, typename P>
Scalar toSeconds(std::chrono::duration<R, P> duration)
{
    using Second = std::chrono::duration<Scalar, std::ratio<1>>;
    return std::chrono::duration_cast<Second>(duration).count();
}

template<typename Scalar, typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
Scalar toSeconds(T duration)
{ return static_cast<Scalar>(duration); }

// overload for nice static_assert message in case of unuspported type
template<typename Scalar, typename T, std::enable_if_t<!std::is_floating_point_v<T>, bool> = true>
Scalar toSeconds(T duration)
{ static_assert(AlwaysFalse<T>::value, "Given type not supported for representation of time values"); }

} // namespace Detail::TimeLoop
#endif // DOXYGEN

/*!
 * \ingroup Core
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
template<class S>
class TimeLoopBase
{
public:
    using Scalar = S;

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
     * \brief Set the current time step size to a given value.
     * \param dt The new value for the time step size \f$\mathrm{[s]}\f$
     */
    template<class Rep, class Period>
    void setTimeStepSize(std::chrono::duration<Rep, Period> dt)
    { setTimeStepSize(Detail::TimeLoop::toSeconds<Scalar>(dt)); }

    /*!
     * \brief Returns true if the simulation is finished.
     */
    virtual bool finished() const = 0;
};

/*!
 * \ingroup Core
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

    template<class Rep1, class Period1,
             class Rep2, class Period2,
             class Rep3, class Period3>
    TimeLoop(std::chrono::duration<Rep1, Period1> startTime,
             std::chrono::duration<Rep2, Period2> dt,
             std::chrono::duration<Rep3, Period3> tEnd,
             bool verbose = true)
    : TimeLoop(
        Detail::TimeLoop::toSeconds<Scalar>(startTime),
        Detail::TimeLoop::toSeconds<Scalar>(dt),
        Detail::TimeLoop::toSeconds<Scalar>(tEnd),
        verbose
    ){}

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
    template<class Rep1, class Period1,
             class Rep2, class Period2,
             class Rep3, class Period3>
    void reset(std::chrono::duration<Rep1, Period1> startTime,
               std::chrono::duration<Rep2, Period2> dt,
               std::chrono::duration<Rep3, Period3> tEnd,
               bool verbose = true)
    {
        reset(
            Detail::TimeLoop::toSeconds<Scalar>(startTime),
            Detail::TimeLoop::toSeconds<Scalar>(dt),
            Detail::TimeLoop::toSeconds<Scalar>(tEnd)
        );
    }

    /*!
     * \brief Reset the time loop
     */
    void reset(Scalar startTime, Scalar dt, Scalar tEnd, bool verbose = true)
    {
        verbose_ =
            verbose &&
            Dune::MPIHelper::getCommunication().rank() == 0;

        startTime_ = startTime;
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
    template<class ScalarOrDuration>
    void setTime(ScalarOrDuration t)
    { time_ = Detail::TimeLoop::toSeconds<Scalar>(t); }

    /*!
     * \brief Set the current simulated time and the time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     * \param stepIdx The new time step index
     */
    template<class ScalarOrDuration>
    void setTime(ScalarOrDuration t, int stepIdx)
    {
        time_ = Detail::TimeLoop::toSeconds<Scalar>(t);
        timeStepIdx_ = stepIdx;
    }

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
            std::cout << Fmt::format("Set new end time to t = {:.5g} seconds.\n", t);
    }

    /*!
     * \brief Returns the current wall clock time (cpu time) spend in this time loop
     */
    double wallClockTime() const
    { return timer_.elapsed(); }

    using TimeLoopBase<Scalar>::setTimeStepSize;
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
        // Warn if dt is so small w.r.t. current time that it renders float addition meaningless
        // For instance, consider (may depend on architecture):
        //     double cien = 100;
        //     double mil = 1000;
        //     if (cien + 1e-14 == cien) std::cout << "Will not be printed" << std::endl;
        //     if (mil + 1e-14 == mil) std::cout << "Will be printed" << std::endl;
        if (!finished() && (time_ + timeStepSize_ == time_))
            std::cerr << Fmt::format("You have set a very small timestep size (dt = {:.5g}).", timeStepSize_)
                      << " This might lead to numerical problems!\n";
    }

    /*!
     * \brief Set the maximum time step size to a given value.
     *
     * \param maxDt The new value for the maximum time step size \f$\mathrm{[s]}\f$
     * \note This may also reduce the currently set timestep size if needed to comply with the set maximum
     */
    template<class ScalarOrDuration>
    void setMaxTimeStepSize(ScalarOrDuration maxDt)
    {
        userSetMaxTimeStepSize_ = Detail::TimeLoop::toSeconds<Scalar>(maxDt);
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
     * \brief The previous time step size
     */
    Scalar previousTimeStepSize() const
    { return previousTimeStepSize_; }

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
        return finished_ || (endTime_ - time_) < baseEps_*(time_ - startTime_);
    }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    {
        return finished() || (endTime_ - time_ - timeStepSize_) < baseEps_*timeStepSize_;
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
            using std::round;
            const auto percent = round( (time_ - startTime_) / (endTime_ - startTime_) * 100 );
            std::cout << Fmt::format("[{:3.0f}%] ", percent)
                      << Fmt::format("Time step {} done in {:.2g} seconds. ", timeStepIdx_, timeStepWallClockTime_)
                      << Fmt::format("Wall clock time: {:.5g}, time: {:.5g}, time step size: {:.5g}\n", cpuTime, time_, previousTimeStepSize_);
        }
    }

    /*!
     * \brief Print final status and stops tracking the time.
     */
    template< class Communicator = Dune::Communication<typename Dune::MPIHelper::MPICommunicator> >
    void finalize(const Communicator& comm = Dune::MPIHelper::getCommunication())
    {
        auto cpuTime = timer_.stop();

        if (verbose_)
            std::cout << Fmt::format("Simulation took {:.5g} seconds on {} processes.\n", cpuTime, comm.size());

        if (comm.size() > 1)
            cpuTime = comm.sum(cpuTime);

        if (verbose_)
            std::cout << Fmt::format("The cumulative CPU time was {:.5g} seconds.\n", cpuTime);
    }

    //! If the time loop has verbose output
    bool verbose() const
    { return verbose_; }

    //! Sets time loop verbosity
    void setVerbose(bool verbose = true)
    { verbose_ = verbose; }

    /*
     * @}
     */

protected:
    static constexpr Scalar baseEps_ = 1e-10;

    Dune::Timer timer_;
    Scalar time_;
    Scalar endTime_;
    Scalar startTime_;

    Scalar timeStepSize_;
    Scalar previousTimeStepSize_;
    Scalar userSetMaxTimeStepSize_;
    Scalar timeAfterLastTimeStep_, timeStepWallClockTime_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};

// always fall back to floating-point representation
template<class Rep1, class Period1,
         class Rep2, class Period2,
         class Rep3, class Period3>
TimeLoop(std::chrono::duration<Rep1, Period1>,
         std::chrono::duration<Rep2, Period2>,
         std::chrono::duration<Rep3, Period3>,
         bool verbose = true) -> TimeLoop<std::conditional_t<
    std::is_floating_point_v<std::common_type_t<Rep1, Rep2, Rep3>>,
    std::common_type_t<Rep1, Rep2, Rep3>,
    double
>>;

/*!
 * \ingroup Core
 * \brief A time loop with a check point mechanism
 */
template <class Scalar>
class CheckPointTimeLoop : public TimeLoop<Scalar>
{
    class CheckPointType {
        static constexpr std::size_t manualIdx = 0;
        static constexpr std::size_t periodicIdx = 1;

    public:
        bool isPeriodic() const { return set_[periodicIdx]; }
        bool isManual() const { return set_[manualIdx]; }
        bool isAny() const { return set_.any(); }

        CheckPointType& withPeriodic(bool value) { set_[periodicIdx] = value; return *this; }
        CheckPointType& withManual(bool value) { set_[manualIdx] = value; return *this; }

    private:
        std::bitset<2> set_;
    };

public:
    template<class... Args>
    CheckPointTimeLoop(Args&&... args)
    : TimeLoop<Scalar>(std::forward<Args>(args)...)
    {
        periodicCheckPoints_ = false;
        deltaPeriodicCheckPoint_ = 0.0;
        lastPeriodicCheckPoint_ = this->startTime_;
        isCheckPoint_ = false;
    }

    /*!
     * \brief Advance time step.
     */
    void advanceTimeStep() override
    {
        const auto dt = this->timeStepSize();
        const auto newTime = this->time()+dt;

        //! Check point management, TimeLoop::isCheckPoint() has to be called after this!
        const auto cpType = nextCheckPointType_(newTime);
        if (cpType.isManual()) checkPoints_.pop();
        if (cpType.isPeriodic()) lastPeriodicCheckPoint_ += deltaPeriodicCheckPoint_;
        isCheckPoint_ = cpType.isAny();

        const auto previousTimeStepSize = this->previousTimeStepSize();

        // advance the time step like in the parent class
        TimeLoop<Scalar>::advanceTimeStep();

        // if this is a check point we might have reduced the time step to reach this check point
        // reset the time step size to the time step size before this time step
        if (!this->willBeFinished())
        {
            using std::max;
            if (isCheckPoint_)
                this->setTimeStepSize(max(dt, previousTimeStepSize));

            // if there is a check point soon check if the time step after the next time step would be smaller
            // than 20% of the next time step, if yes increase the suggested next time step to exactly reach the check point
            // (in the limits of the maximum time step size)
            auto nextDt = this->timeStepSize();
            const auto threshold = 0.2*nextDt;
            const auto nextTime = this->time() + nextDt;

            const auto nextDtToCheckPoint = maxDtToCheckPoint_(nextTime);
            if (nextDtToCheckPoint > Scalar{0} && Dune::FloatCmp::le(nextDtToCheckPoint, threshold))
                nextDt += nextDtToCheckPoint;

            assert(nextDt > 0.0);
            this->setTimeStepSize(nextDt);
        }
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
    template<class ScalarOrDuration1, class ScalarOrDuration2 = Scalar>
    void setPeriodicCheckPoint(ScalarOrDuration1 interval, ScalarOrDuration2 offset = 0.0)
    {
        setPeriodicCheckPoint_(
            Detail::TimeLoop::toSeconds<Scalar>(interval),
            Detail::TimeLoop::toSeconds<Scalar>(offset)
        );
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
    template<class ScalarOrDuration>
    void setCheckPoint(ScalarOrDuration t)
    {
        // set the check point
        setCheckPoint_(Detail::TimeLoop::toSeconds<Scalar>(t));

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

    /*!
     * \brief add checkpoints to the queue from a list of time points
     * \note checkpoints have to be provided in ascending order
     * \param checkPoints the list of check points
     * \note This also updates the time step size and potentially reduces the time step size to meet the next check point
     */
    template<class ScalarOrDuration>
    void setCheckPoint(const std::initializer_list<ScalarOrDuration>& checkPoints)
    { setCheckPoint(checkPoints.begin(), checkPoints.end()); }

    template<class ScalarOrDuration>
    [[deprecated("Use setCheckpoint(begin, end) instead. Will be removed after release 3.9.")]]
    void setCheckPoint(const std::vector<ScalarOrDuration>& checkPoints)
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
            setCheckPoint_(Detail::TimeLoop::toSeconds<Scalar>(*first));

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

private:
    bool fuzzyEqual_(const Scalar t0, const Scalar t1) const
    { return Dune::FloatCmp::eq(t0, t1, this->baseEps_*this->timeStepSize()); }

    void setPeriodicCheckPoint_(Scalar interval, Scalar offset = 0.0)
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
            std::cout << Fmt::format("Enabled periodic check points every {:.5g} seconds ", interval)
                      << Fmt::format("with the next check point at {:.5g} seconds.\n", lastPeriodicCheckPoint_ + interval);

        // check if the current time point is a periodic check point
        if (nextCheckPointType_(this->time() + deltaPeriodicCheckPoint_).isPeriodic())
            isCheckPoint_ = true;

        // make sure we respect this check point on the next time step
        this->setTimeStepSize(this->timeStepSize());
    }

    //! Adds a check point to the queue
    void setCheckPoint_(Scalar t)
    {
        if (Dune::FloatCmp::le(t - this->time(), 0.0, this->timeStepSize()*this->baseEps_))
        {
            if (this->verbose())
                std::cerr << Fmt::format("Couldn't insert checkpoint at t = {:.5g} ", t)
                          << Fmt::format("because that's not in the future! (current simulation time is {:.5g})\n", this->time());
            return;
        }

        if (!checkPoints_.empty())
        {
            if (t <= checkPoints_.back())
            {
                if (this->verbose())
                    std::cerr << Fmt::format("Couldn't insert checkpoint at t = {:.5g} ", t)
                              << Fmt::format("because it's earlier than or equal to the last check point (t = {:.5g}) in the queue.\n", checkPoints_.back())
                              << "Checkpoints can only be inserted in ascending order." << std::endl;
                return;
            }
        }

        checkPoints_.push(t);
        if (this->verbose())
            std::cout << Fmt::format("Set check point at t = {:.5g} seconds.\n", t);
    }

    /*!
     * \brief Return the type of (next) check point at the given time
     */
    CheckPointType nextCheckPointType_(Scalar t)
    {
        return CheckPointType{}
            .withPeriodic(periodicCheckPoints_ && fuzzyEqual_(t - lastPeriodicCheckPoint_, deltaPeriodicCheckPoint_))
            .withManual(!checkPoints_.empty() && fuzzyEqual_(t - checkPoints_.front(), 0.0));
    }

    /*!
     * \brief Aligns dt to the next check point
     */
    Scalar computeStepSizeRespectingCheckPoints_() const
    { return maxDtToCheckPoint_(this->time()); }

    /*!
     * \brief Compute a time step size respecting upcoming checkpoints, starting from the given time t.
     */
    Scalar maxDtToCheckPoint_(Scalar t) const
    {
        static constexpr auto unset = std::numeric_limits<Scalar>::max();
        const auto dtToPeriodic = dtToNextPeriodicCheckPoint_(t);
        const auto dtToManual = dtToNextManualCheckPoint_(t);

        using std::min;
        return min(
            dtToPeriodic.value_or(unset),
            dtToManual.value_or(unset)
        );
    }

    /*!
     * \brief Compute the time step size to the next periodic check point.
     */
    std::optional<Scalar> dtToNextPeriodicCheckPoint_(Scalar t) const
    {
        if (periodicCheckPoints_)
            return lastPeriodicCheckPoint_ + deltaPeriodicCheckPoint_ - t;
        return {};
    }

    /*!
     * \brief Compute a time step size respecting the next manual check point.
     */
    std::optional<Scalar> dtToNextManualCheckPoint_(Scalar t) const
    {
        if (!checkPoints_.empty())
            return checkPoints_.front() - t;
        return {};
    }

    bool periodicCheckPoints_;
    Scalar deltaPeriodicCheckPoint_;
    Scalar lastPeriodicCheckPoint_;
    std::queue<Scalar> checkPoints_;
    bool isCheckPoint_;
};

template<class... Args>
CheckPointTimeLoop(Args&&... args) -> CheckPointTimeLoop<
    typename decltype(TimeLoop{std::forward<Args>(args)...})::Scalar
>;

} // end namespace Dumux

#endif
