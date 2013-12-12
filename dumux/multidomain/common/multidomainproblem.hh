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

#ifndef DUMUX_MULTIDOMAIN_PROBLEM_HH
#define DUMUX_MULTIDOMAIN_PROBLEM_HH

#include "multidomainmodel.hh"
#include "multidomainnewtoncontroller.hh"
#include "multidomainpropertydefaults.hh"
#include "subdomainpropertydefaults.hh"
#include "multidomainassembler.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

/*
* \brief docme
*/

namespace Dumux
{

/*!
 * \ingroup ModelCoupling
 * \brief Base class for problems which involve two sub problems
 *
 * \todo Please docme more!
 */
template<class TypeTag>
class MultiDomainProblem
{
    template<int dim>
    struct VertexLayout
    {
        bool contains(Dune::GeometryType gt) {
            return gt.dim() == 0;
        }
    };

private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, LocalResidual) LocalResidual1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, LocalResidual) LocalResidual2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Problem) SubDomainProblem1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Problem) SubDomainProblem2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, GridView) SubDomainGridView1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, GridView) SubDomainGridView2;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) HostGrid;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MultiDomainGrid;

    typedef typename MultiDomainGrid::LeafGridView MultiDomainGridView;
    typedef typename MultiDomainGrid::Traits::template Codim<0>::Entity MultiDomainElement;
    typedef typename MultiDomainGrid::SubDomainGrid SubDomainGrid;
    typedef typename SubDomainGrid::template Codim<0>::EntityPointer SubDomainElementPointer;

    typedef Dune::MultiDomainMCMGMapper<MultiDomainGridView, VertexLayout> VertexMapper;

public:
    /*!
      * \brief docme
      *
      * \param hostGrid docme
      * \param timeManager The TimeManager which is used by the simulation
      *
      */
    MultiDomainProblem(MultiDomainGrid &mdGrid,
            		   TimeManager &timeManager)
        : timeManager_(timeManager)
		, newtonMethod_(asImp_())
		, newtonCtl_(asImp_())
		, mdGrid_(mdGrid)
		, mdGridView_(mdGrid.leafView())
		, mdVertexMapper_(mdGrid_.leafView())
		, sdID1_(0)
		, sdID2_(1)
		, sdGrid1_(mdGrid.subDomain(sdID1_))
		, sdGrid2_(mdGrid.subDomain(sdID2_))
		, sdProblem1_(timeManager, sdGrid1_.leafView())
		, sdProblem2_(timeManager, sdGrid2_.leafView())
    {  };

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
        sdProblem1().init();
        sdProblem2().init();

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
     * \brief Serialize the simulation's state to disk
     */

    void serialize()
    {
        typedef Dumux::Restart Restarter;
        Restarter res;
        res.serializeBegin(this->asImp_());
        std::cerr << "Serialize to file '" << res.fileName() << "'\n";

        this->timeManager().serialize(res);
        this->asImp_().serialize(res);
        res.serializeEnd();
    }

    /*!
     * \brief Load a previously saved state of the whole simulation
     *        from disk.
     *
     * \param tRestart The simulation time on which the program was
     *                 written to disk.
     */

    void restart(Scalar tRestart)
    {
        typedef Dumux::Restart Restarter;

        Restarter res;

        res.deserializeBegin(this->asImp_(), tRestart);
        std::cout << "Deserialize from file '" << res.fileName() << "'\n";
        this->timeManager().deserialize(res);
        this->asImp_().deserialize(res);
        res.deserializeEnd();
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
        this->model().deserialize(res);
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
        asImp_().sdProblem1().preTimeStep();
        asImp_().sdProblem2().preTimeStep();
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
        asImp_().sdProblem1().postTimeStep();
        asImp_().sdProblem2().postTimeStep();
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
     * very time step. This file is intended to be overwritten by the
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

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {
    	asImp_().sdProblem1().advanceTimeLevel();
    	asImp_().sdProblem2().advanceTimeLevel();

        model_.advanceTimeLevel();
    }

    /*!
     * \brief Write the relevant quantities of the current solution into an VTK output file.
     */
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            asImp_().sdProblem1().writeOutput();
            asImp_().sdProblem2().writeOutput();
        }
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

    /*!
     * \brief Returns the ID of the first domain
     */
    const typename MultiDomainGrid::SubDomainType sdID1() const
    { return sdID1_; }

    /*!
     * \brief Returns the ID of the second domain
     */
    const typename MultiDomainGrid::SubDomainType sdID2() const
    { return sdID2_; }

    // MAY BE THROWN OUT??
    /*!
     * \brief Returns a reference to subproblem1
     */
    SubDomainProblem1& sdProblem1()
    { return sdProblem1_; }

    /*!
     * \brief Returns a reference to subproblem2
     */
    SubDomainProblem2& sdProblem2()
    { return sdProblem2_; }

    /*!
     * \brief Returns a reference to the localresidual1
     */
    LocalResidual1& localResidual1()
    { return sdProblem1().model().localResidual(); };

    /*!
     * \brief Returns a reference to the localresidual2
     */
    LocalResidual2& localResidual2()
    { return sdProblem2().model().localResidual(); };

    /*!
     * \brief Returns a reference to the multidomain grid
     */
    MultiDomainGrid& mdGrid()
    { return mdGrid_; }

    /*!
     * \brief Returns a const reference to the multidomain grid
     */
    const MultiDomainGrid& mdGrid() const
    { return mdGrid_; }

    /*!
     * \brief Returns the multidomain gridview
     */
    const MultiDomainGridView& mdGridView() const
    { return mdGridView_; }


    /*!
     * \brief Returns the multidomain gridview
     */
    const MultiDomainGridView& gridView() const
    { return mdGridView_; }

    /*!
     * \brief Provides a vertex mapper for the multidomain
     *
     */
    VertexMapper& mdVertexMapper()
    { return mdVertexMapper_; }

    /*!
     * \brief Returns a const reference to the subdomain1 grid
     */
    const SubDomainGrid& sdGrid1() const
    { return sdGrid1_; }

    /*!
     * \brief Returns a const reference to the subdomain2 grid
     */
    const SubDomainGrid& sdGrid2() const
    { return sdGrid2_; }

    DUNE_DEPRECATED_MSG("use sdID1 instead")
    const typename MultiDomainGrid::SubDomainType subID1() const
    { return sdID1(); }

    DUNE_DEPRECATED_MSG("use sdID2 instead")
    const typename MultiDomainGrid::SubDomainType subID2() const
    { return sdID2(); }

    DUNE_DEPRECATED_MSG("use sdProblem1 instead")
    SubDomainProblem1& subProblem1()
    { return sdProblem1(); }

    DUNE_DEPRECATED_MSG("use sdProblem2 instead")
    SubDomainProblem2& subProblem2()
    { return sdProblem2(); }

    DUNE_DEPRECATED_MSG("use sdGridView1 instead")
    const SubDomainGridView1 gridView1() const
    { return mdGrid().subDomain(sdID1_).leafView(); }

    DUNE_DEPRECATED_MSG("use sdGridView2 instead")
    const SubDomainGridView2 gridView2() const
    { return mdGrid().subDomain(sdID2_).leafView(); }

    /*!
     * \brief Returns the gridview of subdomain1
     */
    const SubDomainGridView1 sdGridView1() const
    { return sdGrid1_.leafView(); }

    /*!
     * \brief Returns the gridview of subdomain2
     */
    const SubDomainGridView2 sdGridView2() const
    { return sdGrid2_.leafView(); }

    /*!
     * \brief Returns a pointer to the subdomain1 element
     * \param mdElement1 docme
     */
    SubDomainElementPointer sdElementPointer1(const MultiDomainElement& mdElement1)
    { return sdGrid1_.subDomainEntityPointer(mdElement1); }

    /*!
     * \brief Returns a pointer to the subdomain2 element
     *
     * \param mdElement2 docme
     */
    SubDomainElementPointer sdElementPointer2(const MultiDomainElement& mdElement2)
    { return sdGrid2_.subDomainEntityPointer(mdElement2); }


protected:
    /*
    * \brief docme
    * \Returns the implementation of the problem (i.e. static polymorphism)
    */
    void initMortarElements()
    {}


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

	MultiDomainGrid &mdGrid_;
    const MultiDomainGridView mdGridView_;
    VertexMapper mdVertexMapper_;

    typename MultiDomainGrid::SubDomainType sdID1_;
    typename MultiDomainGrid::SubDomainType sdID2_;

    const SubDomainGrid& sdGrid1_;
    const SubDomainGrid& sdGrid2_;

    SubDomainProblem1 sdProblem1_;
    SubDomainProblem2 sdProblem2_;
};

// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag>
std::string MultiDomainProblem<TypeTag>::simname_="simCoupled"; //initialized with default "sim"

}

#endif
