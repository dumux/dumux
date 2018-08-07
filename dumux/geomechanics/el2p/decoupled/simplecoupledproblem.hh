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
#ifndef DUMUX_SIMPLE_COUPLED_PROBLEM_HH
#define DUMUX_SIMPLE_COUPLED_PROBLEM_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/common/timemanager.hh>

namespace Dumux
{

namespace Properties
{
// new type tag for the simple coupling
// NumericModel provides Scalar, GridCreator, ParameterTree
NEW_TYPE_TAG(SimpleCoupled, INHERITS_FROM(NumericModel));

// property forward declarations
NEW_PROP_TAG(GridView);

// property tags that will be set in the problem at hand
NEW_PROP_TAG(SubProblem1TypeTag);
NEW_PROP_TAG(SubProblem2TypeTag);
NEW_PROP_TAG(Problem);

// property tags with default value
NEW_PROP_TAG(TimeManager);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(SolutionVector);

// default property value for the time manager
SET_TYPE_PROP(SimpleCoupled, TimeManager, Dumux::TimeManager<TypeTag>);

// default property value for the grid
SET_PROP(SimpleCoupled, Grid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
public:
    typedef typename GET_PROP_TYPE(SubTypeTag1, Grid) type;
};

// default property value for the solution vector
SET_PROP(SimpleCoupled, SolutionVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
public:
    typedef typename GET_PROP_TYPE(SubTypeTag1, SolutionVector) type;
};

}

/*!
 * \ingroup ModelCoupling
 * \brief Base class for problems which involve two sub problems
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class SimpleCoupledProblem
{
private:
    // the following properties are also required by start.hh, so they are
    // contained in the coupled TypeTag
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // obtain the type tags of the subproblems
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) SubTypeTag2;

    // obtain all other types from the SubTypeTags
    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) SubProblem1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) SubProblem2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, TimeManager) SubTimeManager1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, TimeManager) SubTimeManager2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

public:
    SimpleCoupledProblem(TimeManager &timeManager, const GridView &gridView)
    : timeManager_(timeManager),
    subProblem1_(subTimeManager1_, gridView), subProblem2_(subTimeManager2_, gridView),
    gridView_(gridView)
    {}

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem and the sub-problems.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // set start time for the sub problems
        Scalar tStart = timeManager_.time();

        // set end time for the sub problems
        Scalar tEnd = tStart + timeManager_.timeStepSize();

        bool restart = false;
        // HACK: assume that we restart if time > 0
        if (tStart > 0)
            restart = true;

        // get initial time step size for the subproblems
        Scalar dtSubProblem1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);
        Scalar dtSubProblem2 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);

        // initialize the subproblem time managers
        // (this also initializes the subproblems)
        subTimeManager1_.init(subProblem1_, tStart, dtSubProblem1, tEnd, restart);
        subTimeManager2_.init(subProblem2_, tStart, dtSubProblem2, tEnd, restart);
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
    {}

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        std::cout << "coupled timeIntegration t = " << timeManager_.time() << std::endl;

        // run first model
        subTimeManager1_.setTime(timeManager_.time());
        subTimeManager1_.setEndTime(timeManager_.time() + timeManager_.timeStepSize());
        subTimeManager1_.setTimeStepSize(subTimeManager1_.previousTimeStepSize());
        subTimeManager1_.run();

        // run second model
        subTimeManager2_.setTime(timeManager_.time());
        subTimeManager2_.setEndTime(timeManager_.time() + timeManager_.timeStepSize());
        subTimeManager2_.setTimeStepSize(subTimeManager2_.previousTimeStepSize());
        subTimeManager2_.run();
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {}

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(const Scalar dt)
    {
        return timeManager_.timeStepSize();
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Returns true if the current state of the simulation
     * should be written to disk
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
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It could be either overwritten by the problem files, or simply
     * declared over the setName() function in the application file.
     */
    const char *name() const
    {
        return "simplecoupled";
    }

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {
        subProblem1_.advanceTimeLevel();
        subProblem2_.advanceTimeLevel();
    }

    /*!
     * \brief Write the relevant quantities of the current solution into
     * an VTK output file.
     */
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            subProblem1_.writeOutput();
            subProblem2_.writeOutput();
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
    {}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Returns a reference to subproblem1
     */
    SubProblem1& subProblem1()
    { return subProblem1_; }

    /*!
     * \brief Returns a const reference to subproblem1
     */
    const SubProblem1& subProblem1() const
    { return subProblem1_; }

    /*!
     * \brief Returns a reference to subproblem2
     */
    SubProblem2& subProblem2()
    { return subProblem2_; }

    /*!
     * \brief Returns a const reference to subproblem2
     */
    const SubProblem2& subProblem2() const
    { return subProblem2_; }


    TimeManager &timeManager_;
    SubTimeManager1 subTimeManager1_;
    SubTimeManager2 subTimeManager2_;
private:
    SubProblem1 subProblem1_;
    SubProblem2 subProblem2_;

    const GridView gridView_;
};

}

#endif
