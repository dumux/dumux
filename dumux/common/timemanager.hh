// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Simplify the handling of time dependent problems
 */
#ifndef DUMUX_TIME_MANAGER_HH
#define DUMUX_TIME_MANAGER_HH

#include <dumux/common/propertysystem.hh>

#include <boost/format.hpp>

#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>


namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Problem);
}
/*!
 * \addtogroup SimControl Simulation Supervision
 */

/*!
 * \ingroup SimControl
 * \brief Simplify the handling of time dependent problems.
 *
 * This class manages a sequence of "episodes" which determine the
 * boundary conditions of a problem. This approach is handy if the
 * problem is not static, i.e. that boundary conditions change
 * over time.
 *
 * This class is a low level way to simplify time management for the
 * simulation. It doesn't manage any user data, but only keeps track
 * about what the current "episode" of the simulation is. An episode
 * is a span of simulated time at which the problem behaves in a
 * specific way. It is characerized by the (simulation) time it
 * starts, its length and a consecutive index starting at 0.
 */
template <class TypeTag>
class TimeManager
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
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
        timeStepIdx_ = 0;
        finished_ = false;

        episodeLength_ = 1e100;

        if (verbose_)
            std::cout <<
                "Welcome aboard DuMuX airlines. Please fasten your seatbelts! Emergency exits are near the time integration.\n";
    }

    /*!
     * \brief Initialize the model and problem and write the initial
     *        condition to disk.
     *
     * \param problem The physical problem which needs to be solved
     * \param tStart The start time \f$\mathrm{[s]}\f$ of the simulation (typically 0)
     * \param dtInitial The initial time step size \f$\mathrm{[s]}\f$
     * \param tEnd The time at which the simulation is finished \f$\mathrm{[s]}\f$
     * \param writeInitialSol Specifies whether the initial condition
     *                        should be written to disk
     */
    void init(Problem &problem,
              Scalar tStart,
              Scalar dtInitial,
              Scalar tEnd,
              bool writeInitialSol = true)
    {
        if (verbose_)
            std::cout << "Initializing problem '" << problem_->name() << "'\n";

        problem_ = &problem;
        time_ = tStart;
        timeStepSize_ = dtInitial;
        endTime_ = tEnd;

        // initialize the problem
        problem_->init();

        // initialize the problem
        if (writeInitialSol) {
            time_ -= timeStepSize_;
            problem_->writeOutput();
            time_ += timeStepSize_;
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
    { timeStepSize_ = std::min(dt, maxTimeStepSize()); }

    /*!
     * \brief Returns the suggested time step length \f$\mathrm{[s]}\f$ so that we
     *        don't miss the beginning of the next episode or cross
     *        the end of the simlation.
     */
    Scalar timeStepSize() const
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

        return
            std::min(episodeMaxTimeStepSize(),
                     std::max<Scalar>(0.0, endTime() - time()));
    };

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
     */
    void startNextEpisode(Scalar tStart,
                          Scalar len)
    {
        ++ episodeIndex_;
        episodeStartTime_ = tStart;
        episodeLength_ = len;
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
    bool episodeIsOver() const
    { return time() >= episodeStartTime_ + episodeLength(); }

    /*!
     * \brief Returns true if the current episode will be finished
     *        after the current time step.
     */
    bool episodeWillBeOver() const
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
        if (episodeIsOver())
            return 0.0;

        // make sure that we don't exceed the end of the
        // current episode.
        return
            std::max<Scalar>(0.0,
                             episodeLength() - (time() - episodeStartTime()));
    };

    /*
     * @}
     */

    /*!
     * \brief Runs the simulation using a given problem class.
     *
     * This method makes sure that time steps sizes are aligned to
     * episode boundaries, amongst other stuff.
     */
    void run()
    {
        Dune::Timer timer;
        timer.reset();

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
            if (problem_->shouldWriteOutput() || willBeFinished())
                problem_->writeOutput();

            // perpare the model for the next time integration
            problem_->advanceTimeLevel();

            // advance the simulated time by the current time step size
            time_ += dt;
            ++timeStepIdx_;

            // write restart file if mandated by the problem
            if (problem_->shouldWriteRestartFile())
                problem_->serialize();

            // notify the problem if an episode is finished
            Scalar nextDt = problem_->nextTimeStepSize(dt);
            if (episodeIsOver()) {
                problem_->episodeEnd();
                // the problem sets the initial time step size of a
                // new episode...
                nextDt = timeStepSize();
            }

            // notify the problem that the timestep is done and ask it
            // for a suggestion for the next timestep size
            // set the time step size for the next step
            setTimeStepSize(nextDt);

            if (verbose_) {
                std::cout <<
                    boost::format("Time step %d done. Wall time:%.4g, time:%.4g, time step size:%.4g\n")
                    %timeStepIndex()%timer.elapsed()%time()%dt;
            }

        }

        if (verbose_)
            std::cout << "Simulation took " << timer.elapsed() <<" seconds.\n"
                      << "We hope that you enjoyed simulating with us\n"
                      << "and that you will chose us next time, too.\n";

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
                              << time_ << " "
                              << timeStepIdx_ << " ";
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

    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};
}

#endif
