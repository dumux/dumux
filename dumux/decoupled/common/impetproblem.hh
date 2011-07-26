/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Andreas Lauser                      *
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
#ifndef DUMUX_IMPETPROBLEM_HH
#define DUMUX_IMPETPROBLEM_HH

#include "impetproperties.hh"
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/decoupled/common/gridadapt.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the diffusion problem
 * @author Bernd Flemisch
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::VertexMapper VertexMapper;
    typedef typename SolutionTypes::ElementMapper ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) IMPETModel;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportSolutionType)) TransportSolutionType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureModel)) PressureModel;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportModel)) TransportModel;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
    enum
    {
        wetting = 0, nonwetting = 1,
        adaptiveGrid = GET_PROP_VALUE(TypeTag, PTAG(AdaptiveGrid))
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    // The module to adapt grid. If adaptiveGrid is false, this model does nothing.
    typedef GridAdapt<TypeTag, adaptiveGrid> GridAdaptModel;

    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    //private!! copy constructor
    IMPETProblem(const IMPETProblem&)
    {}

public:
    //! Constructs an object of type IMPETProblemProblem
    /** @param gridView gridview to the grid.
     *  @param verbose verbosity.
     */
    IMPETProblem(const GridView &gridView, bool verbose = true)
    DUNE_DEPRECATED // use IMPETProblem(TimeManager&, const GridView&)
        : gridView_(gridView),
          grid_(0),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          variables_(gridView),
          outputInterval_(1)
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

#if HAVE_MPI
        // communicate to get the bounding box of the whole domain
        for (int i = 0; i < dim; ++i) {
            bboxMin_[i] = gridView.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView.comm().max(bboxMax_[i]);
        };
#endif

        timeManager_ = new TimeManager(verbose);
        deleteTimeManager_ = true;

        pressModel_ = new PressureModel(asImp_());

        transportModel_ = new TransportModel(asImp_());
        model_ = new IMPETModel(asImp_()) ;

        // create an Object to handle adaptive grids
        if (adaptiveGrid)
            gridAdapt_ = new GridAdaptModel(asImp_());

        resultWriter_ = NULL;
    }

    //! Constructs an object of type IMPETProblemProblem
    /** @param gridView gridview to the grid.
     *  @param verbose verbosity.
     */
    IMPETProblem(TimeManager &timeManager, const GridView &gridView)
        : gridView_(gridView),
          grid_(0),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          timeManager_(&timeManager),
          variables_(gridView),
          outputInterval_(1),
          deleteTimeManager_(false)
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

#if HAVE_MPI
        // communicate to get the bounding box of the whole domain
        for (int i = 0; i < dim; ++i) {
            bboxMin_[i] = gridView.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView.comm().max(bboxMax_[i]);
        };
#endif

        pressModel_ = new PressureModel(asImp_());

        transportModel_ = new TransportModel(asImp_());
        model_ = new IMPETModel(asImp_()) ;

        // create an Object to handle adaptive grids
        if (adaptiveGrid)
            gridAdapt_ = new GridAdaptModel(asImp_());

        resultWriter_ = NULL;
    }

    //! destructor
    virtual ~IMPETProblem ()
    {
        delete pressModel_;
        delete transportModel_;
        delete model_;
        delete resultWriter_;
        if (deleteTimeManager_)
            delete timeManager_;
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
        asImp_().generalBoundaryTypes(bcTypes, intersection.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position of the center of the boundary intersection
     */
    void generalBoundaryTypes(BoundaryTypes &bcTypes,
                       const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a generalBoundaryTypes() method.");
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
        asImp_().generalDirichlet(values, intersection.geometry().center());
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
    void generalDirichlet(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a generalDirichlet() method.");
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
        asImp_().generalNeumann(values, intersection.geometry().center());
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
    void generalNeumann(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Neumann conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are neumann, but does not provide "
                   "a neumann() method.");
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
        asImp_().generalSource(values, element.geometry().center());
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
    void generalSource(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    { values = Scalar(0.0);  }

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
        asImp_().generalInitial(values, element.geometry().center());
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
    void generalInitial(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a initial() method.");
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
    {
        // if adaptivity is used, this method adapts the grid.
        // if it is not used, this method does nothing.
        if (adaptiveGrid)
            this->gridAdapt().adaptGrid();
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
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
        typedef TransportSolutionType Solution;
        Solution k1 = asImp_().variables().transportedQuantity();

        Scalar t = timeManager().time();
        Scalar dt = 1e100;

        // obtain the first update and the time step size
        model().update(t, dt, k1);

        //make sure t_old + dt is not larger than tend
        dt = std::min(dt, timeManager().episodeMaxTimeStepSize());

        // check if we are in first TS and an initialDt was assigned
        if (t==0. && timeManager().timeStepSize()!=0.)
        {
#if HAVE_MPI
        dt = this->gridView().comm().min(dt);
#endif

            // check if assigned initialDt is in accordance with dt from first transport step
            if (timeManager().timeStepSize() > dt
#if HAVE_MPI
        && this->gridView().comm().rank() == 0
#endif
            )
                Dune::dwarn << "initial timestep of size " << timeManager().timeStepSize()
                            << "is larger then dt= "<<dt<<" from transport" << std::endl;
            // internally assign next timestep size
            dt = std::min(dt, timeManager().timeStepSize());
        }

        // check maximum allowed time step size
        dt = std::min(dt, asImp_().maxTimeStepSize());

        //make sure the right time-step is used by all processes in the parallel case
#if HAVE_MPI
        dt = this->gridView().comm().min(dt);
#endif

        //assign next tiestep size
        timeManager().setTimeStepSize(dt);

        // explicit Euler: Sat <- Sat + dt*N(Sat)
        asImp_().variables().transportedQuantity() += (k1 *= timeManager().timeStepSize());
    }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * This is used to do some janitorial tasks like writing the
     * current solution to disk.
     */
    void postTimeStep()
    {
        asImp_().pressureModel().updateMaterialLaws();
    };

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
    void setTimeStepSize(Scalar dt)
    { timeManager().setTimeStepSize(dt); }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize(Scalar dt)
    { return timeManager().timeStepSize();}

    /*!
     * \brief Returns the maximum allowed time step size [s]
     *
     * By default this the time step size is unrestricted.
     */
    Scalar maxTimeStepSize() const
    { return std::numeric_limits<Scalar>::infinity(); }

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
            (timeManager().timeStepIndex() % int(5*outputInterval_) == 0);
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
        if (timeManager().timeStepIndex() % outputInterval_ == 0 || timeManager().willBeFinished() || timeManager().episodeWillBeOver())
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
     * the application porblem is declared! If not, the default name "sim"
     * will be used.
     *
     * \param newName The problem's name
     */
    void setName(const char *newName)
    {
        simname_ = newName;
    }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView &gridView() const
    { return gridView_; }

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
        if (!adaptiveGrid)
            Dune::dgrave << "adaptivity module was called despite "
                << "adaptivity is disabled in property system \n;" << adaptiveGrid;

        return *gridAdapt_;
    }

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
    const GlobalPosition &bboxMin() const
    { return bboxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition &bboxMax() const
    { return bboxMax_; }

    //! \name Access functions
    //@{
    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return *timeManager_; }

    //! \copydoc Dumux::IMPETProblem::timeManager()
    const TimeManager &timeManager() const
    { return *timeManager_; }

    /*!
     * \brief Returns variables container
     *
     * This provides access to the important variables that are used in the
     * simulation process, such as pressure, saturation etc.
     */
    Variables& variables ()
    { return variables_; }

    //! \copydoc Dumux::IMPETProblem::variables ()
    const Variables& variables () const
    { return variables_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    IMPETModel &model()
    { return *model_; }

    //! \copydoc Dumux::IMPETProblem::model()
    const IMPETModel &model() const
    { return *model_; }

    /*!
     * \brief Returns the pressure model used for the problem.
     */
    PressureModel &pressureModel()
    { return *pressModel_; }

    //! \copydoc Dumux::IMPETProblem::pressureModel()
    const PressureModel &pressureModel() const
    { return *pressModel_; }

    /*!
     * \brief Returns transport model used for the problem.
     */
    TransportModel &transportModel()
    { return *transportModel_; }

    //! \copydoc Dumux::IMPETProblem::transportModel()
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
     * file.)  See Dumux::Restart for details.
     */
    void serialize()
    {
        typedef Dumux::Restart Restarter;

        Restarter res;
        res.serializeBegin(asImp_());
        std::cerr << "Serialize to file " << res.fileName() << "\n";

        timeManager().serialize(res);
        resultWriter_->serialize(res);
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

        timeManager().deserialize(res);
        resultWriter_->deserialize(res);
        model().deserialize(res);

        res.deserializeEnd();
    };
    // \}

    void addOutputVtkFields()
    {
    }

    //! Write the fields current solution into an VTK output file.
    void writeOutput(bool verbose = true)
    {
        if (verbose && gridView().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        if (!resultWriter_)
            resultWriter_ = new VtkMultiWriter(gridView_, asImp_().name());
        if (GET_PROP_VALUE(TypeTag,PTAG(AdaptiveGrid)))
            resultWriter_->gridChanged();
        resultWriter_->beginWrite(timeManager().time() + timeManager().timeStepSize());
        model().addOutputVtkFields(*resultWriter_);
        asImp_().addOutputVtkFields();
        resultWriter_->endWrite();
    }

    // \}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! Returns the applied VTK-writer for the output
    VtkMultiWriter& resultWriter()
    {
        if (!resultWriter_)
            resultWriter_ = new VtkMultiWriter(gridView_, asImp_().name());
        return *resultWriter_;
    }
    //! \copydoc Dumux::IMPETProblem::resultWriter()
    VtkMultiWriter& resultWriter() const
    {
        if (!resultWriter_)
            resultWriter_ = new VtkMultiWriter(gridView_, asImp_().name());
        return *resultWriter_;
    }

private:
    std::string simname_; // a string for the name of the current simulation,
                                  // which could be set by means of an program argument,
                                 // for example.
    const GridView gridView_;
    // pointer to a possibly adaptive grid.
    Grid* grid_;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    TimeManager *timeManager_;
    bool deleteTimeManager_;

    Variables variables_;

    PressureModel* pressModel_;//!< object including the pressure model
    TransportModel* transportModel_;//!< object including the saturation model
    IMPETModel* model_;

    VtkMultiWriter *resultWriter_;
    int outputInterval_;
    GridAdaptModel* gridAdapt_;
};
}
#endif
