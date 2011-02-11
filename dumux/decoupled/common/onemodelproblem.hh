// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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

#ifndef DUMUX_ONE_MODEL_PROBLEM_HH
#define DUMUX_ONE_MODEL_PROBLEM_HH

#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dumux/common/timemanager.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>


/**
 * @file
 * @brief  Base class for definition of an decoupled diffusion (pressure) or transport problem
 * @author Markus Wolff
 */

namespace Dumux
{

/*! \ingroup IMPET
 *
 * @brief Base class for definition of an decoupled diffusion (pressure) or transport problem
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag, class Implementation>
class OneModelProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef Dumux::TimeManager<TypeTag>      TimeManager;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution Solution;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
    enum
    {
        wetting = 0, nonwetting = 1
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    // private!! copy constructor
    OneModelProblem(const OneModelProblem&)
    {}

public:

    //! Constructs an object of type OneModelProblemProblem
    /*!
     *  \tparam TypeTag The TypeTag
     *  \tparam verbose Output level for Dumux::TimeManager
     */
    OneModelProblem(const GridView &gridView, bool verbose = true)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          timeManager_(verbose),
          variables_(gridView),
          resultWriter_(asImp_().name())
    {
        // calculate the bounding box of the grid view
        VertexIterator vIt = gridView.template begin<dim>();
        const VertexIterator vEndIt = gridView.template end<dim>();
        for (; vIt!=vEndIt; ++vIt) {
            for (int i=0; i<dim; i++) {
                bboxMin_[i] = std::min(bboxMin_[i], vIt->geometry().center()[i]);
                bboxMax_[i] = std::max(bboxMax_[i], vIt->geometry().center()[i]);
            }
        }

        model_ = new Model(asImp_()) ;
    }

    //! destructor
    virtual ~OneModelProblem ()
    {
        delete model_;
    }


    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model().initialize();
    }

    /*!
     * \brief Called by Dumux::TimeManager just before the time
     *        integration.
     */
    void preTimeStep()
    { };

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    { };

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * This is used to do some janitorial tasks like writing the
     * current solution to disk.
     */
    void postTimeStep()
    { };

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {}

    /*!
     * \brief Returns the current time step size [seconds].
     */
    Scalar timeStepSize() const
    { return timeManager_.timeStepSize(); }

    /*!
     * \brief Sets the current time step size [seconds].
     */
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setTimeStepSize(dt); }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(Scalar dt)
    { return timeManager_.timeStepSize();}

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behaviour is to write one restart file every 5 time
     * steps. This file is intented to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        return
            timeManager().timeStepIndex() > 0 &&
            (timeManager().timeStepIndex() % 5 == 0);
    }

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

    void addOutputVtkFields()
    {
        model().addOutputVtkFields(resultWriter_);
    }

    //! Write the fields current solution into an VTK output file.
    void writeOutput()
    {
        if (gridView().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        resultWriter_.beginTimestep(timeManager_.time() + timeManager_.timeStepSize(), gridView());
        asImp_().addOutputVtkFields();
        resultWriter_.endTimestep();
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     */
    void episodeEnd()
    {
        std::cerr << "The end of an episode is reached, but the problem "
                  << "does not override the episodeEnd() method. "
                  << "Doing nothing!\n";
    };

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
     * This function sets the simulation name, which should be called before
     * the application problem is declared! If not, the default name "sim"
     * will be used.
     */
    static void setName(const char *newName)
    {
        simname_ = newName;
    }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView &gridView() const
    { return gridView_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalPosition &bboxMin() const
    { return bboxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition &bboxMax() const
    { return bboxMax_; }

    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return timeManager_; }

    /*!
     * \brief \copybrief Dumux::OneModelProblem::timeManager()
     */
    const TimeManager &timeManager() const
    { return timeManager_; }

    /*!
     * \brief Returns variables object.
     */
    Variables& variables ()
    { return variables_; }

    /*!
     * \brief \copybrief Dumux::OneModelProblem::variables()
     */
    const Variables& variables () const
    { return variables_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model &model()
    { return *model_; }

    /*!
     * \brief \copybrief Dumux::OneModelProblem::model()
     */
    const Model &model() const
    { return *model_; }
    // \}


    /*!
     * \name Restart mechanism
     */
    // \{

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Dumux::Restart for details.
     */
    void serialize()
    {
        typedef Dumux::Restart Restarter;

        Restarter res;
        res.serializeBegin(asImp_());
        std::cerr << "Serialize to file " << res.fileName() << "\n";

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model().serialize(res);

        res.serializeEnd();
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    void deserialize(double t)
    {
        typedef Dumux::Restart Restarter;

        Restarter res;
        res.deserializeBegin(asImp_(), t);
        std::cerr << "Deserialize from file " << res.fileName() << "\n";

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model().deserialize(res);

        res.deserializeEnd();
    };

    // \}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \brief \copybrief Dumux::OneModelProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    VtkMultiWriter& resultWriter()
    {
        return resultWriter_;
    }

    VtkMultiWriter& resultWriter() const
    {
        return resultWriter_;
    }

private:
    static std::string simname_; // a string for the name of the current simulation,
                                  // which could be set by means of an program argument,
                                 // for example.
    const GridView gridView_;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    TimeManager timeManager_;

    Variables variables_;

    Model* model_;

    VtkMultiWriter resultWriter_;
};
// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag, class Implementation>
std::string OneModelProblem<TypeTag, Implementation>::simname_="sim"; //initialized with default "sim"
}
#endif
