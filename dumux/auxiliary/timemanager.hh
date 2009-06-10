/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Simplify the handling of time dependent problems
 */
#ifndef DUNE_TIME_MANAGER_HH
#define DUNE_TIME_MANAGER_HH

#include <boost/format.hpp>

#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>

namespace Dune
{
/*!
 * \addtogroup SimControl Simulation Supervision
 * \brief Simplify the handling of time dependent problems.
 *
 * This class manages a sequence of "episodes" which determine the
 * boundary conditions of a problem. This approach is handy if the
 * problem is not static, i.e. that boundary conditions change
 * over time.
 *
 * This class is a low level way to simplify time management for
 * the simulation. It doesn't manage any user data, but only keeps
 * track about what the current "episode" of the simulation is. An
 * episode is a span of simulated time at which the problem
 * behaves in a specific way. It is characerized by the
 * (simulated) time it starts, its length, an identifier, and a
 * consecutive index starting at 0.
 *
 * \todo Change the time manager to the property system (?) 
 * \todo Remove the episode identifier stuff
 */
template <class EpisodeIdentiferT>
class TimeManager
{
public:
    typedef EpisodeIdentiferT EpisodeIdentifer;

    TimeManager(double endTime = 1e100,
                bool verbose = true)
    {
        wasRestarted_ = false;
        verbose_ =
            verbose &&
            MPIHelper::getCollectiveCommunication().rank() == 0;

        episodeIndex_ = 0;
        episodeStartTime_ = 0;

        time_ = 0.0;
        endTime_ = endTime;
        wasRestarted_ = false;

        stepSize_ = 1.0;
        stepNum_ = 0;
        finished_ = false;

        wasRestarted_ = false;

        episodeLength_ = 1e100;
    }

    /*!
     *  \name Simulated time and time step management
     * @{
     */

    /*!
     * \brief Set the current simulated time, don't change the
     * current time step number.
     */
    void setTime(double t)
    { time_ = t; }

    /*!
     * \brief Set the current simulated time and the time step
     * number.
     */
    void setTime(double t, int stepNum)
    { time_ = t; stepNum_ = stepNum; }

    /*!
     * \brief Let some simulated time pass.
     */
    void proceed(double curDt)
    { time_ += curDt; ++stepNum_; }

    /*!
     * \brief Return the current simulated time.
     */
    double time() const
    { return time_; }

    /*!
     * \brief Returns the number of (simulated) seconds which the simulation runs.
     */
    double endTime() const
    { return endTime_; }

    /*!
     * \brief Set the time of simulated seconds at which the simulation runs.
     */
    void setEndTime(double val) 
    { endTime_ = val; }

    /*!
     * \brief Set the suggested time step size to a fixed value.
     *
     * If the step size would exceed the length of the current
     * episode, the timeStep() method will take care that the
     * step size won't exceed the episode, though.
     */
    void setStepSize(double stepSize)
    {
        stepSize_ = stepSize;
    }

    /*!
     * \brief Returns a suggested timestep length so that we don't
     *        miss the beginning of the next episode.
     */
    double stepSize() const
    {
        // if the current episode is over and the simulation
        // wants to give it some extra time, we will return
        // the time step size it suggested instead of trying
        // to align it to the end of the episode.
        if (episodeIsOver())
            return stepSize_;

        // make sure that we don't exceed the end of the
        // current episode.
        return std::min(episodeLength() - (time() - episodeStartTime()),
                        stepSize_);
    }

    /*!
     * \brief Returns number of time steps which have been
     *        executed since t=0.
     */
    int stepNum() const
    { return stepNum_; }

    /*!
     * \brief Specify whether the simulation is finished
     */
    void setFinished(bool yesno = true)
    { finished_ = yesno; }

    /*!
     * \brief Returns true if the simulation is finished.
     */
    bool finished() const
    { return finished_ || time() >= endTime(); }


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
     * \param tStart Time when the episode began
     * \param len    Length of the episode
     */
    void startNextEpisode(const EpisodeIdentifer &id,
                          double tStart,
                          double len)
    {
        ++ episodeIndex_;
        episodeStartTime_ = tStart;
        episodeLength_ = len;

        episode_ = id;
    }

    /*!
     * \brief Change the current episode of the simulation
     *        assuming that the episode starts at the current
     *        time.
     */
    void startNextEpisode(const EpisodeIdentifer &id,
                          double len = 1e100)
    {
        ++ episodeIndex_;
        episodeStartTime_ = time_;
        episodeLength_ = len;

        episode_ = id;
    }

    /*!
     * \brief Start the next episode, but don't change the episode
     *        identifier.
     *
     * \param len  Length of the episode, infinite if not specified.
     */
    void startNextEpisode(double len = 1e100)
    {
        ++ episodeIndex_;
        episodeStartTime_ = time_;
        episodeLength_ = len;
    }

    /*!
     * \brief Returns the identifier of the current episode
     */
    const EpisodeIdentifer &episode() const
    { return episode_; }

    /*!
     * \brief Returns the identifier of the current episode
     */
    EpisodeIdentifer &episode()
    { return episode_; }

    /*!
     * \brief Returns the index of the current episode.
     *
     * The first episode has the index 0.
     */
    int episodeIndex() const
    { return episodeIndex_; }

    /*!
     * \brief Returns the absolute time when the current episode
     *        started.
     */
    double episodeStartTime() const
    { return episodeStartTime_; }

    /*!
     * \brief Returns the length of the current episode in
     *        simulated time.
     */
    double episodeLength() const
    { return std::min(episodeLength_,  endTime_ - episodeStartTime_); }

    /*!
     * \brief Returns true if the current episode is over.
     */
    bool episodeIsOver() const
    { return time() >= episodeStartTime_ + (1-1e-6)*episodeLength(); }

    /*
     * @}
     */

    /*!
     * \brief Runs the simulation using a given problem class.
     *
     * This method makes sure that time steps sizes are aligned to
     * episode boundaries, amongst other stuff.
     */
    template <class Problem>
    void runSimulation(Problem &problem)
    {
        if (verbose_)
            std::cout <<
                "Welcome aboard DuMuX airlines. Please fasten your seatbelts! Emergency exits are near the time integration.\n";

        Dune::Timer timer;
        timer.reset();

        // initialize them model and write the initial
        // condition to disk
        problem.init();

        // do the time steps
        double curStepSize = 0;
        double nextStepSize = 0;
        while (!finished())
        {
            // execute the time integration (i.e. Runge-Kutta
            // or Euler).
            curStepSize = stepSize();
            problem.timeIntegration(curStepSize, nextStepSize);


            // advance the simulated time by the timestep size
            // actually used by the solver
            proceed(curStepSize);
            // make sure that if don't use an step size
            // aligned to the end of the episode as next step
            // size
            if (!episodeIsOver())
                setStepSize(nextStepSize);

            // notify the problem that the next solution has
            // been found
            problem.timestepDone();

            if (verbose_) {
                std::cout <<
                    boost::format("Timestep %d done: Realtime=%.2f, Simtime=%.4f StepSize=%.4f, NextStepSize=%.2f\n")
                    %stepNum()%timer.elapsed()%time()%curStepSize%stepSize();
            }
        }

        if (verbose_)
            std::cout <<
                boost::format("Simulation took %.3f seconds. Hopefully you enjoyed simulating with us and that you choose us next time as well.\n")
                %timer.elapsed();
    }

    /*!
     * \name Saving/restoring the object state
     * @{
     */
    /*!
     * \brief Write the time manager's state to a restart file.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.serializeSection("TimeManager");
        res.serializeStream() << episodeIndex_ << " "
                              << episodeStartTime_ << " "
                              << time_ << " "
                              << stepNum_ << "\n";
    };

    /*!
     * \brief Read the time manager's state from a restart file.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.deserializeSection("TimeManager");
        res.deserializeStream() >> episodeIndex_
                                >> episodeStartTime_
                                >> time_
                                >> stepNum_;


        std::string dummy;
        std::getline(res.deserializeStream(), dummy);
        wasRestarted_ = true;
    };

    /*
     * @}
     */
private:
    int              episodeIndex_;
    double           episodeStartTime_;
    double           episodeLength_;
    EpisodeIdentifer episode_;

    double time_;
    double endTime_;

    double stepSize_;
    int    stepNum_;
    bool   finished_;
    bool   verbose_;
    bool   wasRestarted_;
};
}

#endif
