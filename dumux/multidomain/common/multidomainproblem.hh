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
 * \brief Base class for problems which involve two sub problems
 */
#ifndef DUMUX_COUPLED_PROBLEM_HH
#define DUMUX_COUPLED_PROBLEM_HH

#include "coupledmodel.hh"
#include "couplednewtoncontroller.hh"

namespace Dumux
{

/*!
 * \ingroup ModelCoupling
 * \brief Base class for problems which involve two sub problems
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class CoupledProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) SubProblem1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) SubProblem2;

public:
    CoupledProblem(TimeManager &timeManager)
        : timeManager_(timeManager),
          newtonMethod_(asImp_()),
          newtonCtl_(asImp_())
    { };

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem and the sub-problems.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // initialize the sub-problems
        asImp_().subProblem1().init();
        asImp_().subProblem2().init();

        // specify the elements which couple
        asImp_().addMetaElements();

        // set the initial condition of the model
        model().init(asImp_());

        // intialize lagrange multipliers
        asImp_().initMortarElements();
    }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * \param res docme
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Dumux::Restart for details.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * \param res docme
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
    }

    /*!
     * \name Simulation control
     */
    // \{

    /*!
     * \brief Start the simulation procedure.
     *
     * \param dtInitial docme
     * \param tEnd docme
     *
     * This method is usually called by the main() function and simply
     * uses Dumux::TimeManager::runSimulation() to do the actual
     * work.
     */
    bool simulate(Scalar dtInitial, Scalar tEnd)
    {
        // set the initial time step and the time where the simulation ends
        timeManager_.setEndTime(tEnd);
        timeManager_.setTimeStepSize(dtInitial);
        timeManager_.runSimulation(asImp_());
        return true;
    };


    /*!
     * \brief Called by the time manager before the time integration. Calls preTimeStep()
     *        of the subproblems.
     */
    void preTimeStep()
    {
        asImp_().subProblem1().preTimeStep();
        asImp_().subProblem2().preTimeStep();
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        // TODO: should be called from the group Implicit
        const int maxFails =
                GET_PARAM_FROM_GROUP(TypeTag, int, Newton, MaxTimeStepDivisions);
        for (int i = 0; i < maxFails; ++i)
        {
            if (model_.update(newtonMethod_, newtonCtl_))
                return;

            // update failed
            Scalar dt = timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            timeManager().setTimeStepSize(nextDt);

            std::cout << "Newton solver did not converge. Retrying with time step of "
                      << timeManager().timeStepSize() << "sec\n";
        }

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " timestep divisions. dt="
                   << timeManager().timeStepSize());
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution. Calls postTimeStep()
     *        of the subproblems.
     */
    void postTimeStep()
    {
        asImp_().subProblem1().postTimeStep();
        asImp_().subProblem2().postTimeStep();
    }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * \param dt docme
     *
     */
    Scalar nextTimeStepSize(const Scalar dt)
    {
        return newtonCtl_.suggestTimeStepSize(dt);
    };

    /*!
     * \brief This method is called by the model if the update to the
     *        next time step failed completely.
     */
    void updateSuccessful()
    {
    	model_.updateSuccessful();
    };

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behaviour is to write out every the solution for
     * very time step. This file is intented to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Returns true if the current state of the simulation should be written to
     *        disk
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \brief Called by the time manager after the end of an episode.
     */
    void episodeEnd()
    {
        std::cerr << "The end of an episode is reached, but the problem "
                  << "does not override the episodeEnd() method. "
                  << "Doing nothing!\n";
    }

    // \}

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It could be either overwritten by the problem files, or simply
     * declared over the setName() function in the application file.
     */
    const char *name() const
    {
        return simname_.c_str();
    }

    /*!
     * \brief Set the problem name.
     *
     * \param newName docme
     *
     * This function sets the simulation name, which should be called before
     * the application porblem is declared! If not, the default name "sim"
     * will be used.
     */
    static void setName(const char *newName)
    {
        simname_ = newName;
    }

    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return timeManager_; }

    /*!
     * \copydoc timeManager()
     */
    const TimeManager &timeManager() const
    { return timeManager_; }

    /*!
     * \brief Returns NewtonControler object used by the simulation
     */
    NewtonController &newtonController()
    { return newtonCtl_; }

    /*!
     * \brief Returns NewtonControler object used by the simulation
     */
    const NewtonController &newtonController() const
    { return newtonCtl_; }

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
    // \}

    void computeCouplingIndices() const
    { DUNE_THROW(Dune::NotImplemented, "Coupled problems need to implement the computeCouplingIndices method!"); }

    SubProblem1 &subProblem1()
    { DUNE_THROW(Dune::NotImplemented, "Coupled problems need to implement the subProblem1 method!"); }
    const SubProblem1 &subProblem1() const
    { DUNE_THROW(Dune::NotImplemented, "Coupled problems need to implement the subProblem1 method!"); }

    SubProblem2 &subProblem2()
    { DUNE_THROW(Dune::NotImplemented, "Coupled problems need to implement the subProblem2 method!"); }
    const SubProblem2 &subProblem2() const
    { DUNE_THROW(Dune::NotImplemented, "Coupled problems need to implement the subProblem2 method!"); }


    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {
    	asImp_().subProblem1().advanceTimeLevel();
    	asImp_().subProblem2().advanceTimeLevel();

        model_.advanceTimeLevel();
    }

    /*!
     * \brief Write the relevant quantities of the current solution into an VTK output file.
     */
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            asImp_().subProblem1().writeOutput();
            asImp_().subProblem2().writeOutput();
        }
    }

    /*!
     * \brief Serialize the simulation's state to disk
     */
    void serialize()
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Dumux::CoupledProblem::serialize");
    }

protected:
    void addMetaElements() 
    {}

    void initMortarElements()
     {}

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    // a string for the name of the current simulation, which could be
    // set by means of an program argument, for example.
    static std::string simname_;

    TimeManager &timeManager_;
    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;
    Model model_;
};
// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag>
std::string CoupledProblem<TypeTag>::simname_="simCoupled"; //initialized with default "sim"

}

#endif
