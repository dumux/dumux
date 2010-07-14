// $Id: timemanager.hh 3736 2010-06-15 09:52:10Z lauser $
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
#ifndef DUMUX_TIME_MANAGER_HH
#define DUMUX_TIME_MANAGER_HH

#include <boost/format.hpp>

#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>

namespace Dumux
{
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
            Dune::MPIHelper::getCollectiveCommunication().rank() == 0;

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
     * \brief Let one time step size pass.
     */
    void proceed()
    { time_ += stepSize_; ++stepNum_; }

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
    void setTimeStepSize(double stepSize)
    {
        stepSize_ = stepSize;
    }

    /*!
     * \brief Returns a suggested timestep length so that we don't
     *        miss the beginning of the next episode.
     */
    double timeStepSize() const
    {
        return stepSize_;
    }

    /*!
     * \brief Returns number of time steps which have been
     *        executed since t=0.
     */
    int timeStepNum() const
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
    { return time() >= episodeStartTime_ + (1 - 1e-14)*episodeLength(); }


    /*!
     * \brief Aligns dt to the episode boundary if t+dt exceeds the current episode.
     */
    double episodeMaxTimeStepSize() const
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
            std::max(0.0,
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
        double dtInitial = stepSize_;
        stepSize_ = 0.0;
        problem.init();
        stepSize_ = dtInitial;

        // do the time steps
        while (!finished())
        {
            problem.timeStepBegin();

            // execute the time integration (i.e. Runge-Kutta
            // or Euler).
            problem.timeIntegration();

            // advance the simulated time by the current time step
            // size
            double dt = timeStepSize();
            proceed();

            problem.timeStepEnd();

            // notify the problem that the timestep is done and ask it
            // for a suggestion for the next timestep size
            double nextDt =
                    std::min(problem.nextTimeStepSize(),
                             episodeMaxTimeStepSize());

            if (verbose_) {
                std::cout <<
                    boost::format("Timestep %d done. CPUt=%.4g, t=%.4g, StepSize=%.4g, NextStepSize=%.4g\n")
                    %timeStepNum()%timer.elapsed()%time()%dt%nextDt;
            }

            // set the time step size for the next step
            setTimeStepSize(nextDt);
        }

        if (verbose_)
            std::cout <<
                boost::format("Simulation took %.3f seconds. Hopefully you enjoyed simulating with us. We hope that you also simulate with us next time.\n")
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
    }

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
    }

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
