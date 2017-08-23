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
 * \ingroup MultiDomain
 * \brief Base class for problems with which involve two sub problems, boundingboxtree-based
 */

#ifndef DUMUX_MULTIDOMAIN_PROBLEM_HH
#define DUMUX_MULTIDOMAIN_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/multidomain/properties.hh>
#include <dumux/multidomain/model.hh>

namespace Dumux
{

/*!
 * \ingroup MultiDomain
 * \brief Base class for probems which involve two sub problems
 */
template<class TypeTag>
class MultiDomainProblem
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, CouplingManager) CouplingManager;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag) StokesProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag) DarcyProblemTypeTag;

    // obtain types from the sub problem type tags
    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, Problem) StokesProblem;
    typedef typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem) DarcyProblem;

    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, TimeManager) StokesTimeManager;
    typedef typename GET_PROP_TYPE(DarcyProblemTypeTag, TimeManager) DarcyTimeManager;

    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, GridView) StokesGridView;
    typedef typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView) DarcyGridView;

    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, SolutionVector) StokesSolutionVector;
    typedef typename GET_PROP_TYPE(DarcyProblemTypeTag, SolutionVector) DarcySolutionVector;

public:
    MultiDomainProblem(TimeManager &timeManager,
                          const StokesGridView &stokesGridView,
                          const DarcyGridView &darcyGridView)
    : timeManager_(timeManager),
      stokesGridView_(stokesGridView),
      darcyGridView_(darcyGridView),
      stokesProblem_(stokesTimeManager_, stokesGridView),
      darcyProblem_(darcyTimeManager_, darcyGridView),
      couplingManager_(std::make_shared<CouplingManager>(stokesProblem_, darcyProblem_)),
      newtonMethod_(asImp_()),
      newtonCtl_(asImp_())
    {
        // check if we got the right dimensions
        static_assert(int(StokesGridView::dimension) == int(DarcyGridView::dimension), "The Stokes grid dimension has to be equal to the Darcy grid dimension.");
        static_assert(int(StokesGridView::dimensionworld) == int(DarcyGridView::dimensionworld), "The subproblems have to live in the same world dimension.");

        // iterative or monolithic
//        useIterativeSolver_ = GET_PARAM_FROM_GROUP(TypeTag, bool, MultiDomain, UseIterativeSolver);

        // read data from the input file
        name_           = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        episodeTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeTime);
        maxTimeStepSize_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, MaxTimeStepSize);

//         if (useIterativeSolver_)
//         {
//             maxIterations_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, IterativeAlgorithm, MaxIterations);
//             tolerance_      = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, IterativeAlgorithm, Tolerance);
//             verbose_    = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, IterativeAlgorithm, Verbose);
//         }
    }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem and the sub-problems.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        Scalar tStart = timeManager().time();
        Scalar tEnd = timeManager().endTime();

        const bool restart = tStart > 0;

        Scalar dtStokesProblem;
        Scalar dtDarcyProblem;

        // get initial time step size for the subproblems
//         if (useIterativeSolver_)
//         {
//             dtStokesProblem = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitialStokesProblem);
//             dtDarcyProblem = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitialDarcyProblem);
//         }
//         else
//         {
            dtStokesProblem = timeManager().timeStepSize();
            dtDarcyProblem = timeManager().timeStepSize();
//         }

        // start new episode
        // episodes are used to keep the subproblems in sync. At the end of each episode
        // of the coupled problem the subproblems need to arrive at the same time.
        // this is relevant if the time step sizes of the subproblems are different.
        timeManager().startNextEpisode(episodeTime_);
        stokesTimeManager().startNextEpisode(episodeTime_);
        darcyTimeManager().startNextEpisode(episodeTime_);

        // initialize the problem coupler prior to the problems
        couplingManager().preInit();

        // set the coupling manager in the sub problems
        stokesProblem().setCouplingManager(couplingManager_);
        darcyProblem().setCouplingManager(couplingManager_);

        // initialize the subproblem time managers (this also initializes the subproblems)
        // the darcy time manager is initialized first as the stokes problem might need some data from the
        // darcy problem for its initialization (e.g. for models using facet coupling with the global caching)
        darcyTimeManager().init(darcyProblem(), tStart, dtDarcyProblem, tEnd, restart);
        stokesTimeManager().init(stokesProblem(), tStart, dtStokesProblem, tEnd, restart);

        // finalize the problem coupler
        couplingManager().postInit();

//         if (!useIterativeSolver_)
            model().init(asImp_());
    }

    /*!
     * \brief This method writes the complete state of the simulation
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Dumux::Restart for details.
     */
    void serialize()
    {}

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void preTimeStep()
    {
//         if (useIterativeSolver_)
//         {
//             previousSolStokes_ = stokesProblem().model().curSol();
//             previousSolDarcy_ = darcyProblem().model().curSol();
//         }
//         else
//         {
            // call low dim problem first because it might grow and the
            // coupling maps might need to be updated before the stokesproblem
            darcyProblem().preTimeStep();
            stokesProblem().preTimeStep();
//         }
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model. Algorithms for the time integration are
     *        implemented here.
     */
    void timeIntegration()
    {
//         if (useIterativeSolver_)
//             asImp_().timeIntegrationIterative();
//         else
            asImp_().timeIntegrationMonolithic();
    }

//     /*!
//      * \brief Time integration for the iterative solver.
//      */
//     void timeIntegrationIterative()
//     {
//         std::cout << std::endl;
//         std::cout << "Start coupled time integration starting at t = " << timeManager().time() << std::endl;
//
//         Scalar error = 1.0;
//         Scalar pStokesNorm = 0.0;
//         Scalar pDarcyNorm = 0.0;
//         StokesSolutionVector solStokesDiff;
//         DarcySolutionVector solDarcyDiff;
//         int iterations = 0;
//
//         if(verbose_)
//             std::cout << "Starting iterative algorithm "
//                       << "with tolerance = " << tolerance_
//                       << " and max iterations = " << maxIterations_ << std::endl;
//
//         // iterative algorithm within each time step
//         while(error > tolerance_)
//         {
//             if(iterations >= maxIterations_)
//                 DUNE_THROW(Dumux::NumericalProblem, "Iterative algorithm didn't converge in " + std::to_string(maxIterations_) + " iterations!");
//
//             // run the stokes model
//             if(verbose_) std::cout << "Solving the stokes problem " << stokesProblem().name() << std::endl;
//             stokesTimeManager().setTime(timeManager().time(), 0);
//             stokesTimeManager().setEndTime(timeManager().time() + timeManager().timeStepSize());
//             stokesTimeManager().setTimeStepSize(stokesTimeManager().previousTimeStepSize());
//             stokesTimeManager().run(); // increases the time step index
//
//             // run the Darcy problem
//             if(verbose_) std::cout << "Solving the low-dimensional problem " << darcyProblem().name() << std::endl;
//             darcyTimeManager().setTime(timeManager().time(), 0);
//             darcyTimeManager().setEndTime(timeManager().time() + timeManager().timeStepSize());
//             darcyTimeManager().setTimeStepSize(darcyTimeManager().previousTimeStepSize());
//             darcyTimeManager().run(); // increases the time step index
//
//             // get current solution
//             currentSolStokes_ = stokesProblem().model().curSol();
//             currentSolDarcy_ = darcyProblem().model().curSol();
//
//             // calculate discrete l2-norm of pressure
//             solStokesDiff = previousSolStokes_;
//             solStokesDiff -= currentSolStokes_;
//             solDarcyDiff = previousSolDarcy_;
//             solDarcyDiff -= currentSolDarcy_;
//
//             pStokesNorm = solStokesDiff.infinity_norm()/currentSolStokes_.infinity_norm();
//             pDarcyNorm = solDarcyDiff.infinity_norm()/currentSolDarcy_.infinity_norm();
//
//             // calculate the error
//             error = pStokesNorm + pDarcyNorm;
//
//             // update the previous solution
//             previousSolStokes_ = currentSolStokes_;
//             previousSolDarcy_ = currentSolDarcy_;
//
//             iterations++;
//
//             if(verbose_)
//             {
//                 std::cout << "Iteration " << iterations << " done "
//                       << "with maximum relative shift in stokes domain = " << pStokesNorm
//                       << " and maximum relative shift in lower-dimensional domain = " << pDarcyNorm
//                       << " | total error = " << error << std::endl;
//             }
//
//         }// end while(error > tolerance_)
//
//         if(verbose_)
//         {
//             std::cout << "Iterative algorithm finished with " << iterations << " iterations "
//                       << "and total error = " << error << std::endl;
//         }
//     }

    /*!
     * \brief Time integration for the monolithic solver.
     */
    void timeIntegrationMonolithic()
    {
        const int maxFails = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, MaxTimeStepDivisions);

        std::cout << std::endl;
        std::cout << "Start coupled time integration starting at t = " << timeManager().time() << std::endl;

        for (int i = 0; i < maxFails; ++i)
        {
            if (model().update(newtonMethod(), newtonController()))
                return;

            const Scalar dt = timeManager().timeStepSize();
            const Scalar nextDt = dt / 2;
            timeManager().setTimeStepSize(nextDt);
            stokesTimeManager().setTimeStepSize(nextDt);
            darcyTimeManager().setTimeStepSize(nextDt);

            // update failed
            std::cout << "Newton solver did not converge with dt="<<dt<<" seconds. Retrying with time step of "
                      << nextDt << " seconds\n";
        }

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " time-step divisions. dt="
                   << timeManager().timeStepSize());
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
//         if (!useIterativeSolver_)
//         {
            // call low dim problem first because it might grow and the
            // coupling maps might need to be updated before the stokesproblem
            darcyProblem().postTimeStep();
            stokesProblem().postTimeStep();
//         }
    }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(const Scalar dt)
    {
        Scalar newDt = newtonCtl_.suggestTimeStepSize(dt);
        stokesTimeManager().setTimeStepSize(newDt);
        darcyTimeManager().setTimeStepSize(newDt);

        return newDt;
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     */
    bool shouldWriteOutput() const
    {
        return (timeManager().episodeWillBeFinished() || timeManager().willBeFinished());
    }

    /*!
     * \brief Returns true if the current state of the simulation
     * should be written to disk
     */
    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Returns the user specified maximum time step size
     *
     * Overload in problem for custom needs.
     */
    Scalar maxTimeStepSize() const
    {
        return maxTimeStepSize_;
    }

    /*!
     * \brief Called by the time manager after the end of an episode.
     */
    void episodeEnd()
    {
        timeManager().startNextEpisode(episodeTime_);
        stokesTimeManager().startNextEpisode(episodeTime_);
        darcyTimeManager().startNextEpisode(episodeTime_);

//         if (!useIterativeSolver_)
//         {
            // set the initial time step size of a an episode to the last real time step size before the episode
            Scalar nextDt = std::max(timeManager().previousTimeStepSize(), timeManager().timeStepSize());
            stokesTimeManager().setTimeStepSize(nextDt);
            darcyTimeManager().setTimeStepSize(nextDt);
//         }
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It could be either overwritten by the problem files, or simply
     * declared over the setName() function in the application file.
     */
    const std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {
        model().advanceTimeLevel();

        asImp_().stokesProblem().advanceTimeLevel();
        asImp_().darcyProblem().advanceTimeLevel();

//         if (!useIterativeSolver_)
//         {
            // foward current time to the subproblems for correct output
            stokesTimeManager().setTime(timeManager().time() + timeManager().timeStepSize(), timeManager().timeStepIndex()+1);
            darcyTimeManager().setTime(timeManager().time() + timeManager().timeStepSize(), timeManager().timeStepIndex()+1);
//         }
    }

    /*!
     * \brief Write the relevant quantities of the current solution into
     * an VTK output file.
     */
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput() && !(timeManager().time() < 0))
        {
            asImp_().stokesProblem().writeOutput();
            asImp_().darcyProblem().writeOutput();
        }
    }

    /*!
     * \brief Load a previously saved state of the whole simulation
     *        from disk.
     *
     * \param tRestart The simulation time on which the program was
     *                 written to disk.
     */
    void restart(const Scalar tRestart)
    {
        DUNE_THROW(Dune::NotImplemented, "Restart multidomain simulations.");
    }

    //! Returns the leafGridView of the Stokes problem
    const StokesGridView &stokesGridView() const
    { return stokesGridView_; }

    //! Returns the leafGridView of the Darcy problem
    const DarcyGridView &darcyGridView() const
    { return darcyGridView_; }

    //! Returns the time manager of the multidomain problem
    const TimeManager &timeManager() const
    { return timeManager_; }

    //! Returns the time manager of the multidomain problem
    TimeManager &timeManager()
    { return timeManager_; }

    //! Returns a reference to the Stokes problem
    StokesProblem &stokesProblem()
    { return stokesProblem_; }

    //! Returns a const reference to the Stokes problem
    const StokesProblem &stokesProblem() const
    { return stokesProblem_; }

    //! Returns a reference to the Darcy problem
    DarcyProblem &darcyProblem()
    { return darcyProblem_; }

    //! Returns a const reference to the Darcy problem
    const DarcyProblem &darcyProblem() const
    { return darcyProblem_; }

    //! Returns a const reference to the Stokes time manager
    const StokesTimeManager &stokesTimeManager() const
    { return stokesTimeManager_; }

    //! Returns a reference to the Stokes time manager
    StokesTimeManager &stokesTimeManager()
    { return stokesTimeManager_; }

    //! Returns a const reference to the Darcy time manager
    const DarcyTimeManager &darcyTimeManager() const
    { return darcyTimeManager_; }

    //! Returns a reference to the Darcy time manager
    DarcyTimeManager &darcyTimeManager()
    { return darcyTimeManager_; }

    //! Return reference to the coupling manager
    CouplingManager &couplingManager()
    { return *couplingManager_; }

    //! Return const reference to the coupling manager
    const CouplingManager &couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model &model()
    { return model_; }

    /*!
     * \copydoc model()
     */
    const Model &model() const
    { return model_; }

    /*!
     * \brief Returns the newton method object
     */
    NewtonMethod &newtonMethod()
    { return newtonMethod_; }

    /*!
     * \copydoc newtonMethod()
     */
    const NewtonMethod &newtonMethod() const
    { return newtonMethod_; }

    /*!
     * \brief Returns the newton contoller object
     */
    NewtonController &newtonController()
    { return newtonCtl_; }

    /*!
     * \copydoc newtonController()
     */
    const NewtonController &newtonController() const
    { return newtonCtl_; }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    TimeManager &timeManager_;

    StokesGridView stokesGridView_;
    DarcyGridView darcyGridView_;

    StokesTimeManager stokesTimeManager_;
    DarcyTimeManager darcyTimeManager_;

    Scalar maxTimeStepSize_;
    StokesProblem stokesProblem_;
    DarcyProblem darcyProblem_;

    StokesSolutionVector previousSolStokes_;
    DarcySolutionVector previousSolDarcy_;

    StokesSolutionVector currentSolStokes_;
    DarcySolutionVector currentSolDarcy_;

    std::shared_ptr<CouplingManager> couplingManager_;

    Model model_;

    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;

    Scalar episodeTime_;

//     bool useIterativeSolver_;

//     // for the iterative algorithm
//     int maxIterations_;
//     Scalar tolerance_;
//     bool verbose_;

    std::string name_;
};

}//end namespace Dumux

#endif
