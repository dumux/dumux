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

#include "model.hh"
#include "newtoncontroller.hh"
#include "propertydefaults.hh"
#include "subdomainpropertydefaults.hh"
#include "assembler.hh"

#include <dumux/io/restart.hh>


namespace Dumux
{

/*!
 * \ingroup ImplicitBaseProblems
 * \ingroup MultidomainModel
 * \brief Base class for problems which involve two sub problems (multidomain problems)s
 */
template<class TypeTag>
class MultiDomainProblem
{

private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, LocalResidual) LocalResidual1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, LocalResidual) LocalResidual2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Problem) SubDomainProblem1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Problem) SubDomainProblem2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, GridView) SubDomainGridView1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, GridView) SubDomainGridView2;

    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MultiDomainGrid;
    typedef typename MultiDomainGrid::LeafGridView MultiDomainGridView;
    typedef typename MultiDomainGrid::Traits::template Codim<0>::Entity MultiDomainElement;
    typedef typename MultiDomainGrid::SubDomainGrid SubDomainGrid;
    typedef typename SubDomainGrid::template Codim<0>::EntityPointer SubDomainElementPointer;

    typedef Dune::MultiDomainMCMGMapper<MultiDomainGridView, Dune::MCMGVertexLayout> VertexMapper;

    // copying a problem is not a good idea
    MultiDomainProblem(const MultiDomainProblem &);

public:
    /*!
     * \brief The problem for the coupling of Stokes and Darcy flow
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    template<class GridView>
    MultiDomainProblem(TimeManager &timeManager,
                       GridView gridView)
        : timeManager_(timeManager)
        , newtonMethod_(asImp_())
        , newtonCtl_(asImp_())
    {
        mdGrid_ = std::make_shared<MultiDomainGrid> (GridCreator::grid());
        mdGridView_ = std::make_shared<MultiDomainGridView> (mdGrid_->leafGridView());
        mdVertexMapper_ = std::make_shared<VertexMapper> (mdGrid_->leafGridView());
        sdProblem1_ = std::make_shared<SubDomainProblem1> (timeManager, mdGrid_->subDomain(sdID1()).leafGridView());
        sdProblem2_ = std::make_shared<SubDomainProblem2> (timeManager, mdGrid_->subDomain(sdID2()).leafGridView());
        maxTimeStepSize_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, MaxTimeStepSize);
    }

    //! \copydoc ImplicitProblem::init()
    void init()
    {
        // initialize the sub-problems
        sdProblem1().init();
        sdProblem2().init();

        // set the initial condition of the model
        model().init(asImp_());

        // initialize Lagrange multipliers
        asImp_().initMortarElements();
    }

    //! \copydoc ImplicitProblem::serialize()
    template <class Restarter>
    void serialize(Restarter &res)
    {
        this->model().serialize(res);
    }

    //! \copydoc ImplicitProblem::serialize()
    void serialize()
    {
        typedef Restart Restarter;
        Restarter res;
        res.serializeBegin(this->asImp_());
        std::cout << "Serialize to file '" << res.fileName() << "'\n";
        this->timeManager().serialize(res);
        this->asImp_().serialize(res);
        res.serializeEnd();
    }

    //! \copydoc ImplicitProblem::restart()
    void restart(Scalar tRestart)
    {
        typedef Restart Restarter;
        Restarter res;
        res.deserializeBegin(this->asImp_(), tRestart);
        std::cout << "Deserialize from file '" << res.fileName() << "'\n";
        this->timeManager().deserialize(res);
        this->asImp_().deserialize(res);
        res.deserializeEnd();
    }

    //! \copydoc ImplicitProblem::deserialize()
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
     * \brief Called by the time manager before the time integration. Calls preTimeStep()
     *        of the subproblems.
     */
    void preTimeStep()
    {
        asImp_().sdProblem1().preTimeStep();
        asImp_().sdProblem2().preTimeStep();
    }

    //! \copydoc ImplicitProblem::timeIntegration()
    void timeIntegration()
    {
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

    //! \copydoc ImplicitProblem::maxTimeStepSize()
    Scalar maxTimeStepSize() const
    {
        return maxTimeStepSize_;
    }

    //! \copydoc ImplicitProblem::nextTimeStepSize()
    Scalar nextTimeStepSize(const Scalar dt)
    {
        return newtonCtl_.suggestTimeStepSize(dt);
    }

    /*!
     * \brief This method is called by the model if the update to the
     *        next time step failed completely.
     */
    void updateSuccessful()
    {
        model_.updateSuccessful();
    }

    //! \copydoc ImplicitProblem::shouldWriteOutput()
    bool shouldWriteOutput() const
    { return true; }

    //! \copydoc ImplicitProblem::shouldWriteRestartFile()
    bool shouldWriteRestartFile() const
    { return false; }

    //! \copydoc ImplicitProblem::episodeEnd()
    void episodeEnd()
    {
        std::cerr << "The end of an episode is reached, but the problem "
                  << "does not override the episodeEnd() method. "
                  << "Doing nothing!\n";
    }

    //! \copydoc ImplicitProblem::advanceTimeLevel()
    void advanceTimeLevel()
    {
        asImp_().sdProblem1().advanceTimeLevel();
        asImp_().sdProblem2().advanceTimeLevel();

        model_.advanceTimeLevel();
    }

    //! \copydoc ImplicitProblem::writeOutput()
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            asImp_().sdProblem1().writeOutput();
            asImp_().sdProblem2().writeOutput();
        }
    }


    // \}

    //! \copydoc ImplicitProblem::name()
    const char *name() const
    { return simname_.c_str(); }

    //! \copydoc ImplicitProblem::setName()
    static void setName(const char *newName)
    { simname_ = newName; }

    //! \copydoc ImplicitProblem::timeManager()
    TimeManager &timeManager()
    { return timeManager_; }

    //! \copydoc ImplicitProblem::timeManager()
    const TimeManager &timeManager() const
    { return timeManager_; }

    //! \copydoc ImplicitProblem::newtonController()
    NewtonController &newtonController()
    { return newtonCtl_; }

    //! \copydoc ImplicitProblem::newtonController()
    const NewtonController &newtonController() const
    { return newtonCtl_; }

    //! \copydoc ImplicitProblem::model()
    Model &model()
    { return model_; }

    //! \copydoc ImplicitProblem::model()
    const Model &model() const
    { return model_; }
    // \}

    /*!
     * \brief Returns the ID of the first domain
     */
    const typename MultiDomainGrid::SubDomainIndex sdID1() const
    { return typename MultiDomainGrid::SubDomainIndex(0); }

    /*!
     * \brief Returns the ID of the second domain
     */
    const typename MultiDomainGrid::SubDomainIndex sdID2() const
    { return typename MultiDomainGrid::SubDomainIndex(1); }

    /*!
     * \brief Returns a reference to subproblem1
     */
    SubDomainProblem1& sdProblem1()
    { return *sdProblem1_; }

    /*!
     * \brief Returns a const reference to subproblem1
     */
    const SubDomainProblem1& sdProblem1() const
    { return *sdProblem1_; }

    /*!
     * \brief Returns a reference to subproblem2
     */
    SubDomainProblem2& sdProblem2()
    { return *sdProblem2_; }

    /*!
     * \brief Returns a const reference to subproblem2
     */
    const SubDomainProblem2& sdProblem2() const
    { return *sdProblem2_; }

    /*!
     * \brief Returns a reference to the localresidual1
     */
    LocalResidual1& localResidual1()
    { return sdProblem1().model().localResidual(); }

    /*!
     * \brief Returns a reference to the localresidual2
     */
    LocalResidual2& localResidual2()
    { return sdProblem2().model().localResidual(); }

    /*!
     * \brief Returns a reference to the multidomain grid
     */
    MultiDomainGrid& mdGrid()
    { return *mdGrid_; }

    /*!
     * \brief Returns a const reference to the multidomain grid
     */
    const MultiDomainGrid& mdGrid() const
    { return *mdGrid_; }

    /*!
     * \brief Returns the multidomain gridview
     */
    const MultiDomainGridView& mdGridView() const
    { return *mdGridView_; }

    /*!
     * \brief Returns the multidomain gridview
     */
    const MultiDomainGridView& gridView() const
    { return *mdGridView_; }

    /*!
     * \brief Provides a vertex mapper for the multidomain
     */
    VertexMapper& mdVertexMapper()
    { return *mdVertexMapper_; }

    /*!
     * \brief Returns a const reference to the subdomain1 grid
     */
    const SubDomainGrid& sdGrid1() const
    { return mdGrid_->subDomain(sdID1()); }

    /*!
     * \brief Returns a const reference to the subdomain2 grid
     */
    const SubDomainGrid& sdGrid2() const
    { return mdGrid_->subDomain(sdID2()); }

    /*!
     * \brief Returns the gridview of subdomain1
     */
    const SubDomainGridView1 sdGridView1() const
    { return sdGrid1().leafGridView(); }

    /*!
     * \brief Returns the gridview of subdomain2
     */
    const SubDomainGridView2 sdGridView2() const
    { return sdGrid2().leafGridView(); }

    /*!
     * \brief Returns a pointer to the subdomain1 element
     *
     * \param mdElement1 The multi domain element1
     */
    SubDomainElementPointer sdElementPointer1(const MultiDomainElement& mdElement1)
    { return sdGrid1().subDomainEntityPointer(mdElement1); }

    /*!
     * \brief Returns a pointer to the subdomain2 element
     *
     * \param mdElement2 The multi domain element2
     */
    SubDomainElementPointer sdElementPointer2(const MultiDomainElement& mdElement2)
    { return sdGrid2().subDomainEntityPointer(mdElement2); }


protected:
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
    Scalar maxTimeStepSize_;
    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;

    Model model_;

    std::shared_ptr<MultiDomainGrid> mdGrid_;
    std::shared_ptr<MultiDomainGridView> mdGridView_;
    std::shared_ptr<VertexMapper> mdVertexMapper_;

    std::shared_ptr<SubDomainProblem1> sdProblem1_;
    std::shared_ptr<SubDomainProblem2> sdProblem2_;
};

// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag>
std::string MultiDomainProblem<TypeTag>::simname_ = "simCoupled";

} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_PROBLEM_HH
