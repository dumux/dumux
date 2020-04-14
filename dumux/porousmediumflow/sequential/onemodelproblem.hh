// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

#ifndef DUMUX_ONE_MODEL_PROBLEM_HH
#define DUMUX_ONE_MODEL_PROBLEM_HH

#include <dune/common/shared_ptr.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>


/**
 * @file
 * @brief  Base class for definition of an sequential diffusion (pressure) or transport problem
 */

namespace Dumux
{

/*! \ingroup IMPET
 *
 * @brief Base class for definition of an sequential diffusion (pressure) or transport problem
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class OneModelProblem
{
private:
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;

    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

    using VtkMultiWriter = Dumux::VtkMultiWriter<GridView>;

    using Variables = GetPropType<TypeTag, Properties::Variables>;

    using Model = GetPropType<TypeTag, Properties::Model>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using VertexMapper = typename SolutionTypes::VertexMapper;
    using ElementMapper = typename SolutionTypes::ElementMapper;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum
    {
        wetting = 0, nonwetting = 1
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    // private!! copy constructor
    OneModelProblem(const OneModelProblem&)
    {}

public:

    //! Constructs an object of type OneModelProblemProblem
    /*!
     *  \tparam TypeTag The TypeTag
     *  \tparam verbose Output level for TimeManager
     */
    OneModelProblem(Grid& grid, bool verbose = true)
        : gridView_(grid.leafGridView()),
          bBoxMin_(std::numeric_limits<double>::max()),
          bBoxMax_(-std::numeric_limits<double>::max()),
          variables_(grid.leafGridView()),
          outputInterval_(1),
          outputTimeInterval_(0)
    {
        // calculate the bounding box of the grid view
        using std::max;
        using std::min;
        for (const auto& vertex : vertices(grid.leafGridView())) {
            for (int i=0; i<dim; i++) {
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().center()[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().center()[i]);
            }
        }

        timeManager_ = std::make_shared<TimeManager>(verbose);

        model_ = std::make_shared<Model>(asImp_()) ;
        maxTimeStepSize_ = getParam<Scalar>("TimeManager.MaxTimeStepSize", std::numeric_limits<Scalar>::max());
    }

    //! Constructs an object of type OneModelProblemProblem
    /*!
     *  \tparam TypeTag The TypeTag
     *  \tparam verbose Output level for TimeManager
     */
    OneModelProblem(TimeManager& timeManager, Grid& grid)
        : gridView_(grid.leafGridView()),
          bBoxMin_(std::numeric_limits<double>::max()),
          bBoxMax_(-std::numeric_limits<double>::max()),
          variables_(grid.leafGridView()),
          outputInterval_(1),
          outputTimeInterval_(0)
    {
        // calculate the bounding box of the grid view
        using std::max;
        using std::min;
        for (const auto& vertex : vertices(grid.leafGridView())) {
            for (int i=0; i<dim; i++) {
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().center()[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().center()[i]);
            }
        }

        timeManager_ = Dune::stackobject_to_shared_ptr<TimeManager>(timeManager);

        model_ = std::make_shared<Model>(asImp_()) ;
        maxTimeStepSize_ = getParam<Scalar>("TimeManager.MaxTimeStepSize", std::numeric_limits<Scalar>::max());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param intersection The intersection for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &bcTypes,
                       const Intersection &intersection) const
    {
        // forward it to the method which only takes the global coordinate
        asImp_().boundaryTypesAtPos(bcTypes, intersection.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position of the center of the boundary intersection
     */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes,
                       const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a boundaryTypesAtPos() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param intersection The boundary intersection
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values,
            const Intersection &intersection) const
    {
        // forward it to the method which only takes the global coordinate
        asImp_().dirichletAtPos(values, intersection.geometry().center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position of the center of the boundary intersection
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a dirichletAtPos() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param intersection The boundary intersection
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
            const Intersection &intersection) const
    {
        // forward it to the interface with only the global position
        asImp_().neumannAtPos(values, intersection.geometry().center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param globalPos The position of the center of the boundary intersection
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Neumann conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are neumann, but does not provide "
                   "a neumannAtPos() method.");
    }

    /*!
     * \brief Evaluate the source term
     *
     * \param values The source and sink values for the conservation equations
     * \param element The element
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element) const
    {
        // forward to generic interface
        asImp_().sourceAtPos(values, element.geometry().center());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void sourceAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {         // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a sourceAtPos() method.");
        }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The element
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element) const
    {
        // forward to generic interface
        asImp_().initialAtPos(values, element.geometry().center());
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // initialize with 0 by default
        values = 0;
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        variables_.initialize();
        model().initialize();
    }

    /*!
     * \brief Called by TimeManager just before the time
     *        integration.
     */
    void preTimeStep()
    {}

    /*!
     * \brief Called by TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {}

    /*!
     * \brief Called by TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * This is used to do some janitorial tasks like writing the
     * current solution to disk.
     */
    void postTimeStep()
    {}

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {}

     /*!
     * \brief Returns the user specified maximum time step size
     *
     * Overload in problem for custom needs.
     */
    Scalar maxTimeStepSize() const
    { return maxTimeStepSize_; }

    /*!
     * \brief Returns the current time step size [seconds].
     */
    Scalar timeStepSize() const
    { return timeManager().timeStepSize(); }

    /*!
     * \brief Sets the current time step size [seconds].
     */
    void setTimeStepSize(Scalar dt)
    { timeManager().setTimeStepSize(dt); }

    /*!
     * \brief Called by TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(Scalar dt)
    { return timeManager().timeStepSize();}

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
     * \brief Sets a time interval for Output
     *
     * The default is 0.0 -> Output determined by output number interval (<tt>setOutputInterval(int)</tt>)
     */
    void setOutputTimeInterval(const Scalar timeInterval)
    {
        outputTimeInterval_ = (timeInterval > 0.0) ? timeInterval : 1e100;
        timeManager().startNextEpisode(outputTimeInterval_);
    }

    /*!
     * \brief Sets the interval for Output
     *
     * The default is 1 -> Output every time step
     */
    void setOutputInterval(int interval)
    { outputInterval_ = interval; }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behaviour is to write out every the solution for
     * very time step. This file is intented to be overwritten by the
     * implementation.
     */

    bool shouldWriteOutput() const
    {
        if (outputInterval_ > 0)
        {
            if (timeManager().timeStepIndex() % outputInterval_ == 0
                || timeManager().willBeFinished()
                || timeManager().episodeWillBeFinished())
            {
                return true;
            }
        }
        else if (timeManager().willBeFinished()
                 || timeManager().episodeWillBeFinished() || timeManager().timeStepIndex() == 0)
        {
            return true;
        }
        return false;
    }

    void addOutputVtkFields()
    {}

    //! Write the fields current solution into an VTK output file.
    void writeOutput(bool verbose = true)
    {
        if (verbose && gridView().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView(), asImp_().name());
        resultWriter_->beginWrite(timeManager().time() + timeManager().timeStepSize());
        model().addOutputVtkFields(*resultWriter_);
        asImp_().addOutputVtkFields();
        resultWriter_->endWrite();
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     */
    void episodeEnd()
    {
        if (outputTimeInterval_ > 0.0 && !timeManager().finished())
        {
            timeManager().startNextEpisode(outputTimeInterval_);
        }
        else if (!timeManager().finished())
        {
            std::cerr << "The end of an episode is reached, but the problem "
                      << "does not override the episodeEnd() method. "
                      << "Doing nothing!\n";
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
    const std::string& name() const
    {
        return simname_;
    }

    /*!
     * \brief Set the problem name.
     *
     * This function sets the simulation name, which should be called before
     * the application problem is declared! If not, the default name "sim"
     * will be used.
     */
    void setName(const std::string& newName)
    {
        simname_ = newName;
    }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView &gridView() const
    { return gridView_; }

    /*!
     * \brief Returns the mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return variables_.vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return variables_.elementMapper(); }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return *timeManager_; }

    /*!
     * \brief \copybrief OneModelProblem::timeManager()
     */
    const TimeManager &timeManager() const
    { return *timeManager_; }

    /*!
     * \brief Returns variables object.
     */
    Variables& variables ()
    { return variables_; }

    /*!
     * \brief \copybrief OneModelProblem::variables()
     */
    const Variables& variables () const
    { return variables_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model &model()
    { return *model_; }

    /*!
     * \brief \copybrief OneModelProblem::model()
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
     * file.)  See Restart for details.
     */
    void serialize()
    {
        using Restarter = Restart;

        Restarter res;
        res.serializeBegin(asImp_());
        std::cout << "Serialize to file " << res.fileName() << "\n";

        timeManager().serialize(res);
        resultWriter().serialize(res);
        res.template deserializeEntities<0> (model(), gridView_);

        res.serializeEnd();
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    void restart(double tRestart)
    {
        using Restarter = Restart;

        Restarter res;
        res.deserializeBegin(asImp_(), tRestart);
        std::cout << "Deserialize from file " << res.fileName() << "\n";

        timeManager().deserialize(res);
        resultWriter().deserialize(res);
        res.template deserializeEntities<0> (model(), gridView_);

        res.deserializeEnd();
    }

    // \}

protected:
    VtkMultiWriter& resultWriter()
    {
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView_, asImp_().name());
        return *resultWriter_;
    }

    VtkMultiWriter& resultWriter() const
    {
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView_, asImp_().name());
        return *resultWriter_;
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \brief \copybrief OneModelProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::string simname_; // a string for the name of the current simulation,
                                  // which could be set by means of an program argument,
                                 // for example.
    const GridView gridView_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<TimeManager> timeManager_;
    Scalar maxTimeStepSize_;

    Variables variables_;

    std::shared_ptr<Model> model_;

    std::shared_ptr<VtkMultiWriter> resultWriter_;
    int outputInterval_;
    Scalar outputTimeInterval_;
};

}
#endif
