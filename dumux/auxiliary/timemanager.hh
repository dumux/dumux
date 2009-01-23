/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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

namespace Dune
{
    /*!
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
     */
    template <class EpisodeIdentiferT>
    class TimeManager
    {
    public:
        typedef EpisodeIdentiferT EpisodeIdentifer;

        TimeManager(const EpisodeIdentifer &id,
                    double len = 1e100,
                    bool verbose = true)
            {
                verbose_ = verbose;

                init_();

                episode_ = id;
                episodeLength_ = len;
            }

        TimeManager(double len,
                    bool verbose = true)
            {
                verbose_ = verbose;

                init_();
                episodeLength_ = len;
            }

        TimeManager(bool verbose = true)
            {
                verbose_ = verbose;
                init_();
                episodeLength_ = 1e100;
            }

        /*!
         *  \defgroup time Simulated time and time step management
         * @{
         */

        /*!
         * \brief Set the current simulated time, don't change the
         * current time step number.
         */
        void setTime(double t)
            { currentTime_ = t; }

        /*!
         * \brief Set the current simulated time and the time step
         * number.
         */
        void setTime(double t, int stepNum)
            { currentTime_ = t; stepNum_ = stepNum; }

        /*!
         * \brief Let some simulated time pass.
         */
        void proceed(double curDt)
            { currentTime_ += curDt; ++stepNum_; }

        /*!
         * \brief Return the current simulated time.
         */
        double time() const
            { return currentTime_; }

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
            { return finished_; }


        /*!
         * @}
         */


        /*!
         * \defgroup episode Episode management
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
                episodeStartTime_ = currentTime_;
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
                episodeStartTime_ = currentTime_;
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
            { return episodeLength_; }

        /*!
         * \brief Returns true if the current episode is over.
         */
        bool episodeIsOver() const
            { return time() + episodeLength_*1e-6 >= episodeStartTime_ + episodeLength_; }

        /*!
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
                typedef typename Problem::DomainTraits::Scalar  Scalar;

                Dune::Timer timer;
                timer.reset();


                // reset the time manager
                init_();

                // initialize them model and write the initial
                // condition to disk
                problem.init();

                // do the time steps
                Scalar curStepSize = 0;
                Scalar nextStepSize = 0;
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
                            boost::format("Timestep %d done: Realtime=%.2f, Simtime=%.2f StepSize=%.2g, NextStepSize=%.2g\n")
                            %stepNum()%timer.elapsed()%time()%curStepSize%stepSize();
                    }
                }

                if (verbose_)
                    std::cout <<
                        boost::format("Simulation took %.3f seconds\n")
                        %timer.elapsed();
            }


    private:
        void init_()
            {
                episodeIndex_ = 0;
                episodeStartTime_ = 0;
                currentTime_ = 0.0;
                stepSize_ = 1.0;
                stepNum_ = 0;
                finished_ = false;
            }

        int              episodeIndex_;
        double           episodeStartTime_;
        double           episodeLength_;
        EpisodeIdentifer episode_;

        double currentTime_;
        double stepSize_;
        int    stepNum_;
        bool   finished_;
        bool   verbose_;
    };
}

#endif
