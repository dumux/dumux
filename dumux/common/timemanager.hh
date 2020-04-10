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
#ifndef DUMUX_TIME_MANAGER_HH
#define DUMUX_TIME_MANAGER_HH

#warning "This file is deprecated. Use the TimeLoop class in common/timeloop.hh"

#include <algorithm>

#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>

#include "properties.hh"
#include "parameters.hh"

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
 */
template <class TypeTag>
class TimeManager
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    TimeManager(const TimeManager&)
    {}
public:

    TimeManager(bool verbose = true)
    {
        verbose_ =
            verbose &&
            Dune::MPIHelper::getCollectiveCommunication().rank() == 0;

        episodeIndex_ = 0;
        episodeStartTime_ = 0;

        time_ = 0.0;
        endTime_ = -1e100;

        timeStepSize_ = 1.0;
        previousTimeStepSize_ = timeStepSize_;
        timeStepIdx_ = 0;
        finished_ = false;

        episodeLength_ = 1e100;
    }

    /*!
     * \brief Initialize the model and problem and write the initial
     *        condition to disk.
     *
     * \param problem The physical problem which needs to be solved
     * \param tStart The start time \f$\mathrm{[s]}\f$ of the simulation (typically 0)
     * \param dtInitial The initial time step size \f$\mathrm{[s]}\f$
     * \param tEnd The time at which the simulation is finished \f$\mathrm{[s]}\f$
     * \param restart Specifies whether the initial condition should be written to disk
     */
    void init(Problem &problem,
              Scalar tStart,
              Scalar dtInitial,
              Scalar tEnd,
              bool restart = false)
    {
        timer_.reset();

        problem_ = &problem;
        time_ = tStart;
        timeStepSize_ = dtInitial;
        previousTimeStepSize_ = dtInitial;
        endTime_ = tEnd;

        if (verbose_)
            std::cout << "Initializing problem '" << problem_->name() << "'\n";

        // initialize the problem
        problem_->init();

        // restart problem if necessary
        if(restart)
            problem_->restart(tStart);
        else
        {
            // write initial condition (if problem is not restarted)
            time_ -= timeStepSize_;
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();
            time_ += timeStepSize_;
        }

        if (verbose_) {
            int numProcesses = Dune::MPIHelper::getCollectiveCommunication().size();
            std::cout << "Initialization took " << timer_.elapsed() <<" seconds on "
                      << numProcesses << " processes.\n"
                      << "The cumulative CPU time was " << timer_.elapsed()*numProcesses << " seconds.\n";
        }

    }

    /*!
     *  \name Simulated time and time step management
     * @{
     */

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
    {  return timer_.elapsed(); }

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
     * \param yesno If true the simulation is considered finished
     *              before the end time is reached, else it is only
     *              considered finished if the end time is reached.
     */
    void setFinished(bool yesno = true)
    { finished_ = yesno; }

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
        return min(min(episodeMaxTimeStepSize(),
                       problem_->maxTimeStepSize()),
                   max<Scalar>(0.0, endTime() - time()));
    }

    /*
     * @}
     */


    /*!
     * \name episode Episode management
     * @{
     */

    /*!
     * \brief Change the current episode of the simulation.
     *
     * \param tStart Time when the episode began \f$\mathrm{[s]}\f$
     * \param len    Length of the episode \f$\mathrm{[s]}\f$
     * \param description descriptive name of the episode
     */
    void startNextEpisode(Scalar tStart,
                          Scalar len,
                          const std::string &description = "")
    {
        ++ episodeIndex_;
        episodeStartTime_ = tStart;
        episodeLength_ = len;
        episodeDescription_ = description;
    }

    /*!
     * \brief Start the next episode, but don't change the episode
     *        identifier.
     *
     * \param len  Length of the episode \f$\mathrm{[s]}\f$, infinite if not specified.
     */
    void startNextEpisode(Scalar len = 1e100)
    {
        ++ episodeIndex_;
        episodeStartTime_ = time_;
        episodeLength_ = len;
    }

    /*!
     * \brief Returns the index of the current episode.
     *
     * The first episode has the index 0.
     */
    int episodeIndex() const
    { return episodeIndex_; }

    /*!
     * \brief Returns the absolute time when the current episode
     *        started \f$\mathrm{[s]}\f$.
     */
    Scalar episodeStartTime() const
    { return episodeStartTime_; }

    /*!
     * \brief Returns the length of the current episode in
     *        simulated time \f$\mathrm{[s]}\f$.
     */
    Scalar episodeLength() const
    { return episodeLength_; }

    /*!
     * \brief Returns true if the current episode is finished at the
     *        current time.
     */
    bool episodeIsFinished() const
    { return time() >= episodeStartTime_ + episodeLength(); }

    /*!
     * \brief Returns true if the current episode will be finished
     *        after the current time step.
     */
    bool episodeWillBeFinished() const
    { return time() + timeStepSize() >= episodeStartTime_ + episodeLength(); }

    /*!
     * \brief Aligns the time step size to the episode boundary if the
     *        current time step crosses the boundary of the current episode.
     */
    Scalar episodeMaxTimeStepSize() const
    {
        // if the current episode is over and the simulation
        // wants to give it some extra time, we will return
        // the time step size it suggested instead of trying
        // to align it to the end of the episode.
        if (episodeIsFinished())
            return 0.0;

        // make sure that we don't exceed the end of the
        // current episode.
        using std::max;
        return max<Scalar>(0.0, episodeLength() - (time() - episodeStartTime()));
    }

    /*
     * @}
     */

    /*!
     * \brief Runs the simulation using a given problem class.
     *
     * This method makes sure that time step sizes are aligned to
     * episode boundaries, amongst other stuff.
     */
    void run()
    {
        timer_.reset();

        // do the time steps
        while (!finished())
        {
            // pre-process the current solution
            problem_->preTimeStep();

            // execute the time integration scheme
            problem_->timeIntegration();
            Scalar dt = timeStepSize();

            // post-process the current solution
            problem_->postTimeStep();

            // write the result to disk
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();

            // prepare the model for the next time integration
            problem_->advanceTimeLevel();

            // write restart file if mandated by the problem
            if (problem_->shouldWriteRestartFile())
                problem_->serialize();

            // advance the simulated time by the current time step size
            time_ += dt;
            ++timeStepIdx_;

            if (verbose_) {
                std::cout
                    << "Time step "<<timeStepIndex()<<" done. "
                    << "Wall time:"<<timer_.elapsed()
                    <<", time:"<<time()
                    <<", time step size:"<<dt
                    <<"\n";
            }

            // notify the problem if an episode is finished
            if (episodeIsFinished()) {
                //define what to do at the end of an episode in the problem
                problem_->episodeEnd();

                //check if a time step size was explicitly defined in problem->episodeEnd()
                if (Dune::FloatCmp::eq<Scalar>(dt, timeStepSize()))
                {
                    // set the initial time step size of a an episode to the last real time step size before the episode
                    using std::max;
                    Scalar nextDt = max(previousTimeStepSize_, timeStepSize());
                    previousTimeStepSize_ = nextDt;
                    setTimeStepSize(nextDt);
                }
            }
            else
            {
                // notify the problem that the time step is done and ask it
                // for a suggestion for the next time step size
                // set the time step size for the next step
                previousTimeStepSize_ = timeStepSize();
                setTimeStepSize(problem_->nextTimeStepSize(dt));
            }
        }

        if (verbose_) {
            int numProcesses = Dune::MPIHelper::getCollectiveCommunication().size();
            std::cout << "Simulation took " << timer_.elapsed() <<" seconds on "
                      << numProcesses << " processes.\n"
                      << "The cumulative CPU time was " << timer_.elapsed()*numProcesses << " seconds.\n";
        }
    }

    /*!
     * \name Saving/restoring the object state
     * @{
     */
    /*!
     * \brief Write the time manager's state to a restart file.
     *
     * \tparam Restarter The type of the object which takes care to serialize data
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.serializeSectionBegin("TimeManager");
        res.serializeStream() << episodeIndex_ << " "
                              << episodeStartTime_ << " "
                              << episodeLength_ << " "
                              << time_ + timeStepSize() << " "
                              << timeStepIdx_ + 1 << " ";
        res.serializeSectionEnd();
    }

    /*!
     * \brief Read the time manager's state from a restart file.
     *
     * \tparam Restarter The type of the object which takes care to deserialize data
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.deserializeSectionBegin("TimeManager");
        res.deserializeStream() >> episodeIndex_
                                >> episodeStartTime_
                                >> episodeLength_
                                >> time_
                                >> timeStepIdx_;
        res.deserializeSectionEnd();
    }

    /*
     * @}
     */

private:
    Problem *problem_;
    int episodeIndex_;
    Scalar episodeStartTime_;
    Scalar episodeLength_;
    std::string episodeDescription_;

    Dune::Timer timer_;
    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    Scalar previousTimeStepSize_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};
}

#endif
