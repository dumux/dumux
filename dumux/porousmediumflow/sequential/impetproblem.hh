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
#ifndef DUMUX_IMPETPROBLEM_HH
#define DUMUX_IMPETPROBLEM_HH

#include <dune/common/float_cmp.hh>
#include "impetproperties.hh"
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>
#include <dumux/io/adaptivegridrestart.hh>

#include <dumux/porousmediumflow/sequential/gridadapt.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the diffusion problem
 */

namespace Dumux
{
/*!
 * \ingroup IMPET
 * @brief base class for problems using a sequential implicit-explicit strategy
 *
 *  \tparam TypeTag      problem TypeTag
 *  \tparam Implementation problem implementation
 */
template<class TypeTag>
class IMPETProblem
{
private:
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

    using VtkMultiWriter = Dumux::VtkMultiWriter<GridView>;

    using Variables = GetPropType<TypeTag, Properties::Variables>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using VertexMapper = typename SolutionTypes::VertexMapper;
    using ElementMapper = typename SolutionTypes::ElementMapper;

    using IMPETModel = GetPropType<TypeTag, Properties::Model>;
    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;
    using PressureModel = GetPropType<TypeTag, Properties::PressureModel>;
    using TransportModel = GetPropType<TypeTag, Properties::TransportModel>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    enum
        {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld
        };
    enum
        {
            wetting = 0, nonwetting = 1,
            adaptiveGrid = getPropValue<TypeTag, Properties::AdaptiveGrid>()
        };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    // The module to adapt grid. If adaptiveGrid is false, this model does nothing.
    using GridAdaptModel = GetPropType<TypeTag, Properties::GridAdaptModel>;

    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    //private!! copy constructor
    IMPETProblem(const IMPETProblem&)
    {}

public:
    /*! \brief Constructs an object of type IMPETProblemProblem
     * \param timeManager The time manager
     * \param grid The grid
     */
    IMPETProblem(TimeManager &timeManager, Grid& grid)
        : IMPETProblem(timeManager, grid, grid.leafGridView())
    {}

    IMPETProblem(TimeManager &timeManager, Grid& grid, const GridView& gridView)
        : grid_(&grid),
          gridView_(gridView),
          bBoxMin_(std::numeric_limits<double>::max()),
          bBoxMax_(-std::numeric_limits<double>::max()),
          timeManager_(&timeManager),
          outputInterval_(1),
          outputTimeInterval_(0.0),
          vtkOutputLevel_(-1)
    {
        // calculate the bounding box of the grid view
        using std::min;
        using std::max;
        for (const auto& vertex : vertices(this->gridView())) {
            for (int i=0; i<dim; i++) {
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().center()[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().center()[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        if (this->gridView().comm().size() > 1)
            for (int i = 0; i < dim; ++i) {
                bBoxMin_[i] = this->gridView().comm().min(bBoxMin_[i]);
                bBoxMax_[i] = this->gridView().comm().max(bBoxMax_[i]);
            }

        variables_ = std::make_shared<Variables>(this->gridView());
        pressModel_ = std::make_shared<PressureModel>(asImp_());

        transportModel_ = std::make_shared<TransportModel>(asImp_());
        model_ = std::make_shared<IMPETModel>(asImp_()) ;

        // create an Object to handle adaptive grids
        if (adaptiveGrid)
            gridAdapt_ = std::make_shared<GridAdaptModel>(asImp_());

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
        dtVariationRestrictionFactor_ = getParam<Scalar>("Impet.DtVariationRestrictionFactor", std::numeric_limits<Scalar>::max());
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
        variables_->initialize();
        model().initialize();
        if (adaptiveGrid)
            gridAdapt().init();
    }

    /*!
     * \brief Called by TimeManager just before the time
     *        integration.
     */
    void preTimeStep()
    {
        // if adaptivity is used, this method adapts the grid.
        // if it is not used, this method does nothing.
        if (adaptiveGrid && timeManager().timeStepIndex() > 0)
            this->gridAdapt().adaptGrid();
        asImp_().pressureModel().updateMaterialLaws();
    }

    /*!
     * \brief Called by TimeManager in order to do a time
     *        integration on the model.
     *
     * \note \a timeStepSize and \a nextStepSize are references and may
     *       be modified by the timeIntegration(). On exit of this
     *       function \a timeStepSize must contain the step size
     *       actually used by the time integration for the current
     *       steo, and \a nextStepSize must contain a suggestion for the
     *       next time step size.
     */
    void timeIntegration()
    {
        // allocate temporary vectors for the updates
        TransportSolutionType updateVector;//(this->transportModel().solutionSize());

        Scalar t = timeManager().time();
        Scalar dt = std::numeric_limits<Scalar>::max();

        // obtain the first update and the time step size
        model().update(t, dt, updateVector);

        using std::max;
        using std::min;
        if (t > 0 && timeManager().timeStepSize()> 0.)
        {
            Scalar oldDt = timeManager().timeStepSize();
            Scalar minDt = max(oldDt - dtVariationRestrictionFactor_ * oldDt, 0.0);
            Scalar maxDt = oldDt + dtVariationRestrictionFactor_ * oldDt;
            dt = max(min(dt, maxDt), minDt);
        }

        //make sure t_old + dt is not larger than tend
        dt = min(dt, timeManager().episodeMaxTimeStepSize());

        // check if we are in first TS and an initialDt was assigned
        if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(t, 0.0, 1.0e-30)
            && Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(timeManager().timeStepSize(), 0.0, 1.0e-30))
        {
            if (gridView().comm().size() > 1)
                dt = gridView().comm().min(dt);

            // check if assigned initialDt is in accordance with dt from first transport step
            if (timeManager().timeStepSize() > dt
                && gridView().comm().rank() == 0)
                Dune::dwarn << "Initial timestep of size " << timeManager().timeStepSize()
                            << " is larger then dt=" << dt <<" from transport" << std::endl;
            // internally assign next timestep size
            dt = min(dt, timeManager().timeStepSize());
        }

        //make sure the right time-step is used by all processes in the parallel case
        if (gridView().comm().size() > 1)
            dt = gridView().comm().min(dt);

        // check maximum allowed time step size
        timeManager().setTimeStepSize(dt);

        // explicit Euler: Sat <- Sat + dt*N(Sat)
        transportModel().updateTransportedQuantity(updateVector);
    }

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
     * \brief Returns the current time step size [seconds].
     */
    Scalar timeStepSize() const
    { return timeManager().timeStepSize(); }

    /*!
     * \brief Sets the current time step size [seconds].
     */
    void setTimeStepSize(const Scalar dt)
    { timeManager().setTimeStepSize(dt); }

    /*!
     * \brief Called by TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(const Scalar dt)
    { return timeManager().timeStepSize();}

    /*!
     * \brief Returns the maximum allowed time step size [s]
     *
     * By default this the time step size is unrestricted.
     */
    Scalar maxTimeStepSize() const
    { return maxTimeStepSize_; }

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
        if (outputInterval_ > 0)
        {
            return
                timeManager().timeStepIndex() > 0 &&
                (timeManager().timeStepIndex() % int(100*outputInterval_) == 0);
        }
        else
            return
                shouldWriteOutput();
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
    void setOutputInterval(const int interval)
    {
        using std::max;
        outputInterval_ = max(interval, 0);
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
     * the application porblem is declared! If not, the default name "sim"
     * will be used.
     *
     * \param newName The problem's name
     */
    void setName(std::string newName)
    {
        simname_ = newName;
    }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView& gridView() const
    {
        return gridView_;
    }

    /*!
     * \brief Returns the current grid which used by the problem.
     */
    Grid &grid()
    {
        if (!grid_)
        {
            DUNE_THROW(Dune::InvalidStateException, "Grid was called in problemclass, "
                       << "although it is not specified. Do so by using setGrid() method!");
        }
        return *grid_;
    }
    /*!
     * \brief Specifies the grid from outside the problem.
     * \param grid The grid used by the problem.
     */
    void setGrid(Grid &grid)
    {
        grid_ = &grid;
    }

    /*!
     * \brief Returns adaptivity model used for the problem.
     */
    GridAdaptModel& gridAdapt()
    {
        return *gridAdapt_;
    }

    /*!
     * \brief Capability to introduce problem-specific routines at the
     * beginning of the grid adaptation
     *
     * Function is called at the beginning of the standard grid
     * modification routine, GridAdapt::adaptGrid() .
     */
    void preAdapt()
    {
    }

    /*!
     * \brief Capability to introduce problem-specific routines after grid adaptation
     *
     * Function is called at the end of the standard grid
     * modification routine, GridAdapt::adaptGrid() , to allow
     * for problem-specific output etc.
     */
    void postAdapt()
    {
        // write out new grid
        //        Dune::VTKWriter<LeafGridView> vtkwriter(problem_.gridView());
        //        vtkwriter.write("latestgrid",Dune::VTK::appendedraw);
    }

    /*!
     * \brief Returns the mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return variables_->vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return variables_->elementMapper(); }

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

    //! \name Access functions
    //@{
    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return *timeManager_; }

    //! \copydoc IMPETProblem::timeManager()
    const TimeManager &timeManager() const
    { return *timeManager_; }

    /*!
     * \brief Returns variables container
     *
     * This provides access to the important variables that are used in the
     * simulation process, such as pressure, saturation etc.
     */
    Variables& variables ()
    { return *variables_; }

    //! \copydoc IMPETProblem::variables ()
    const Variables& variables () const
    { return *variables_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    IMPETModel &model()
    { return *model_; }

    //! \copydoc IMPETProblem::model()
    const IMPETModel &model() const
    { return *model_; }

    /*!
     * \brief Returns the pressure model used for the problem.
     */
    PressureModel &pressureModel()
    { return *pressModel_; }

    //! \copydoc IMPETProblem::pressureModel()
    const PressureModel &pressureModel() const
    { return *pressModel_; }

    /*!
     * \brief Returns transport model used for the problem.
     */
    TransportModel &transportModel()
    { return *transportModel_; }

    //! \copydoc IMPETProblem::transportModel()
    const TransportModel &transportModel() const
    { return *transportModel_; }
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

        // do the actual serialization process: write primary variables
        res.template serializeEntities<0> (*pressModel_, gridView());
        res.template serializeEntities<0> (*transportModel_, gridView());

        res.serializeEnd();

        if (adaptiveGrid)
        {
            AdaptiveGridRestart<Grid>::serializeGrid(asImp_());
        }
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     * @param tRestart Restart time
     */
    void restart(const double tRestart)
    {
        if (adaptiveGrid)
        {
            AdaptiveGridRestart<Grid>::restartGrid(asImp_());
            variables().initialize();
            model().initialize();
        }

        using Restarter = Restart;

        Restarter res;
        res.deserializeBegin(asImp_(), tRestart);
        std::cout << "Deserialize from file " << res.fileName() << "\n";

        timeManager().deserialize(res);

        resultWriter().deserialize(res);

        // do the actual serialization process: get primary variables
        res.template deserializeEntities<0> (*pressModel_, gridView());
        res.template deserializeEntities<0> (*transportModel_, gridView());
        pressureModel().updateMaterialLaws();

        res.deserializeEnd();
    }

    void addOutputVtkFields()
    {
    }
    /*!
     * \brief Returns the vtk output verbosity level
     *
     * Level is set by property or input file.
     */
    int vtkOutputLevel() const
    {
        return vtkOutputLevel_;
    }

    //! Write the fields current solution into an VTK output file.
    void writeOutput(bool verbose = true)
    {
        if (verbose && gridView().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView(), asImp_().name());
        if (adaptiveGrid)
            resultWriter_->gridChanged();
        resultWriter_->beginWrite(timeManager().time() + timeManager().timeStepSize());
        model().addOutputVtkFields(*resultWriter_);
        asImp_().addOutputVtkFields();
        resultWriter_->endWrite();
    }

    // \}

protected:
    //! Returns the applied VTK-writer for the output
    VtkMultiWriter& resultWriter()
    {
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView(), asImp_().name());
        return *resultWriter_;
    }
    //! \copydoc IMPETProblem::resultWriter()
    VtkMultiWriter& resultWriter() const
    {
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView(), asImp_().name());
        return *resultWriter_;
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::string simname_; // a string for the name of the current simulation,
    // which could be set by means of a program argument, for example.

    // pointer to a possibly adaptive grid.
    Grid *grid_;
    GridView gridView_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    TimeManager *timeManager_;
    Scalar maxTimeStepSize_;

    std::shared_ptr<Variables> variables_;

    std::shared_ptr<PressureModel> pressModel_;//!< object including the pressure model
    std::shared_ptr<TransportModel> transportModel_;//!< object including the saturation model
    std::shared_ptr<IMPETModel> model_;

    std::shared_ptr<VtkMultiWriter> resultWriter_;
    int outputInterval_;
    Scalar outputTimeInterval_;
    int vtkOutputLevel_;
    std::shared_ptr<GridAdaptModel> gridAdapt_;

    Scalar dtVariationRestrictionFactor_;
};
}
#endif
