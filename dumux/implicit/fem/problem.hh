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
 * \brief Base class for all fully implicit problems
 */
#ifndef DUMUX_IMPLICIT_FEM_PROBLEM_HH
#define DUMUX_IMPLICIT_FEM_PROBLEM_HH

#include "properties.hh"
#include "model.hh"

#include <dumux/io/restart.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitBaseProblems
 * \brief Base class for finite-element fully implicit problems
 *
 * \note All quantities are specified assuming a threedimensional
 *       world. Problems discretized using 2D grids are assumed to be
 *       extruded by \f$1 m\f$ and 1D grids are assumed to have a
 *       cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class ImplicitFemProblem
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using NewtonMethod = typename GET_PROP_TYPE(TypeTag, NewtonMethod);
    using VtkMultiWriter = typename GET_PROP_TYPE(TypeTag, VtkMultiWriter);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! copying a problem is not a good idea
    ImplicitFemProblem(const ImplicitFemProblem&) = delete;

    /*!
     * \brief Constructor
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    ImplicitFemProblem(TimeManager &timeManager, const GridView &gridView)
    : gridView_(gridView),
      bBoxMin_(std::numeric_limits<double>::max()),
      bBoxMax_(-std::numeric_limits<double>::max()),
      elementMapper_(gridView),
      vertexMapper_(gridView),
      timeManager_(&timeManager),
      newtonMethod_(asImp_()),
      newtonCtl_(asImp_()),
      gravity_(0)
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : vertices(gridView))
        {
            for (int i = 0; i < dimWorld; i++)
            {
                bBoxMin_[i] = std::min(bBoxMin_[i], vertex.geometry().corner(0)[i]);
                bBoxMax_[i] = std::max(bBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        if (gridView.comm().size() > 1)
        {
            for (int i = 0; i < dimWorld; ++i)
            {
                bBoxMin_[i] = gridView.comm().min(bBoxMin_[i]);
                bBoxMax_[i] = gridView.comm().max(bBoxMax_[i]);
            }
        }

        // set a default name for the problem
        simName_ = "sim";
        maxTimeStepSize_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, MaxTimeStepSize);

        // if gravity is enabled, set gravitational acceleration
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            gravity_[dimWorld-1]  = -9.81;

        // instantiate the spatial parameters
        spatialParams_ = std::make_shared<SpatialParams>(asImp_(), gridView);
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // set the initial condition of the model
        model().init(asImp_());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param intersection The intersection which is to be checked
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const Intersection& intersection) const
    {
        // forward it to the method which only takes the global coordinate
        return asImp_().boundaryTypesAtPos(intersection.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position of the finite volume in global coordinates
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any boundaryTypes method
        //! set Dirichlet boundary conditions everywhere for all primary variables
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param element The finite element
     * \param intersection The intersection which is to be checked
     * \param ipGlobal the global coordinate of the integration point
     *
     * The method returns the boundary types information.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const Intersection& intersection,
                               const GlobalPosition& ipGlobal) const
    {
        // forward it to the method which only takes the global coordinate
        return asImp_().dirichletAtPos(ipGlobal);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a dirichlet() or a dirichletAtPos() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param intersection The intersection between element and boundary
     * \param elemSol The solution at the dofs of the element
     * \param ipData Data on the shape values and gradients at the integration point
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    PrimaryVariables neumann(const Element& element,
                             const Intersection& intersection,
                             const ElementSolutionVector& elemSol,
                             const IpData& ipData) const
    {
        // forward it to the interface with only the global position
        return asImp_().neumannAtPos(ipData.ipGlobal());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the boundary face's integration point in global coordinates
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any neumann method
        //! return no-flow Neumann boundary conditions at all Neumann boundaries
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param ipData Data on the shape values and gradients at the integration point
     * \param secVars The primary/secondary variables evaluated at the integration point
     *
     */
    PrimaryVariables source(const Element& element,
                            const IpData& ipData,
                            const SecondaryVariables& secVars) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(ipData.ipGlobal());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    PrimaryVariables sourceAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any source method
        //! return 0.0 (no source terms)
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "an initial() or an initialAtPos() method.");
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor(const Element& element,
                           const PrimaryVariables& priVars) const
    {
        // forward to generic interface
        return asImp_().extrusionFactorAtPos(element.geometry().center());
    }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any extrusion factor method
        //! return 1.0
        return 1.0;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     */
    const GlobalPosition &gravityAtPos(const GlobalPosition &pos) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>ProblemEnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void preTimeStep()
    {}

    /*!
     * \brief Called by TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        const int maxFails =
                GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, MaxTimeStepDivisions);
        for (int i = 0; i < maxFails; ++i) {
            if (model_.update(newtonMethod_, newtonCtl_))
                return;

            Scalar dt = timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            timeManager().setTimeStepSize(nextDt);

            // update failed
            std::cout << "Newton solver did not converge with dt="<<dt<<" seconds. Retrying with time step of "
                      << nextDt << " seconds\n";
        }

        // if the simulation  run is about to abort, write restart files for the current and previous time steps:
        // write restart file for the current time step
        serialize();

        //write restart file for the previous time step:
        //set the time manager and the solution vector to the previous time step
        const Scalar time = timeManager().time();
        timeManager().setTime(time - timeManager().previousTimeStepSize());
        const auto curSol = model_.curSol();
        model_.curSol() = model_.prevSol();
        //write restart file
        serialize();
        //reset time manager and solution vector
        model_.curSol() = curSol;
        timeManager().setTime(time);

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " time-step divisions. dt="
                   << timeManager().timeStepSize()
                   << ".\nThe solutions of the current and the previous time steps "
                   << "have been saved to restart files.");
    }

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

    /*!
     * \brief Called by TimeManager whenever a solution for a
     *        time step has been computed and the simulation time has
     *        been updated.
     *
     * \param dt The current time-step size
     */
    Scalar nextTimeStepSize(const Scalar dt)
    { return newtonCtl_.suggestTimeStepSize(dt); }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 5 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    { return timeManager().timeStepIndex() > 0 && (timeManager().timeStepIndex() % 10 == 0); }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for
     * every time step. This function is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    { return true; }

     /*!
     * \brief Returns the user specified maximum time step size
     *
     * Overload in problem for custom needs.
     */
    Scalar maxTimeStepSize() const
    { return maxTimeStepSize_; }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {}

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    {
        model_.advanceTimeLevel();
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
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
    const std::string& name() const
    { return simName_; }

    /*!
     * \brief Set the problem name.
     *
     * This static method sets the simulation name, which should be
     * called before the application problem is declared! If not, the
     * default name "sim" will be used.
     *
     * \param newName The problem's name
     */
    void setName(const std::string& newName)
    {
        simName_ = newName;
    }


    /*!
     * \brief Returns the number of the current VTK file.
     */
    int currentVTKFileNumber()
    {
        createResultWriter_();
        return resultWriter_->curWriterNum();
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
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \brief Returns the mapper for vertices to indices for constant grids.
     */
    const VertexMapper &vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const ElementMapper &elementMapper() const
    { return elementMapper_; }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    VertexMapper &vertexMapper()
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for possibly adaptive grids.
     */
    ElementMapper &elementMapper()
    { return elementMapper_; }

    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return *timeManager_; }

    /*!
     * \copydoc timeManager()
     */
    const TimeManager &timeManager() const
    { return *timeManager_; }

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
     * \name Restart mechanism
     */
    // \{

    /*!
     * \brief This method writes the complete state of the simulation
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Restart for details.
     */
    void serialize()
    {
        typedef Restart Restarter;
        Restarter res;
        res.serializeBegin(asImp_());
        if (gridView().comm().rank() == 0)
            std::cout << "Serialize to file '" << res.fileName() << "'\n";

        timeManager().serialize(res);
        asImp_().serialize(res);
        res.serializeEnd();
    }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Restart for details.
     *
     * \tparam Restarter The serializer type
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        createResultWriter_();
        resultWriter_->serialize(res);
        model().serialize(res);
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
        typedef Restart Restarter;

        Restarter res;

        res.deserializeBegin(asImp_(), tRestart);
        if (gridView().comm().rank() == 0)
            std::cout << "Deserialize from file '" << res.fileName() << "'\n";
        timeManager().deserialize(res);
        asImp_().deserialize(res);
        res.deserializeEnd();
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        createResultWriter_();
        resultWriter_->deserialize(res);
        model().deserialize(res);
    }

    // \}

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by writeOutput().
     */
    void addOutputVtkFields()
    {}

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into an VTK output file.
     */
    void writeOutput(const bool verbose = true)
    {
        // write the current result to disk
        if (verbose && gridView().comm().rank() == 0)
            std::cout << "Writing result file for \"" << asImp_().name() << "\"\n";

        // calculate the time _after_ the time was updated
        Scalar t = timeManager().time() + timeManager().timeStepSize();
        createResultWriter_();
        resultWriter_->beginWrite(t);
        model().addOutputVtkFields(model().curSol(), *resultWriter_);
        asImp_().addOutputVtkFields();
        resultWriter_->endWrite();
    }

    /*!
     * \brief Returns a reference to the grid
     */
    Grid &grid()
    { return GridCreator::grid(); }

    /*!
     * \brief Capability to introduce problem-specific routines at the
     * beginning of the grid adaptation
     *
     * Function is called at the beginning of the standard grid
     * modification routine, GridAdapt::adaptGrid() .
     */
    void preAdapt()
    {}

    /*!
     * \brief Capability to introduce problem-specific routines after grid adaptation
     *
     * Function is called at the end of the standard grid
     * modification routine, GridAdapt::adaptGrid() , to allow
     * for problem-specific output etc.
     */
    void postAdapt()
    {}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! Returns the applied VTK-writer for the output
    VtkMultiWriter& resultWriter()
    {
        createResultWriter_();
        return *resultWriter_;
    }
    //! \copydoc IMPETProblem::resultWriter()
    VtkMultiWriter& resultWriter() const
    {
        createResultWriter_();
        return *resultWriter_;
    }


private:
    // makes sure that the result writer exists
    void createResultWriter_()
    {
        if (!resultWriter_)
            resultWriter_ = std::make_shared<VtkMultiWriter>(gridView_, asImp_().name());
    }

    std::string simName_;
    const GridView gridView_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    TimeManager *timeManager_;
    Scalar maxTimeStepSize_;

    Model model_;

    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;

    GlobalPosition gravity_;
    std::shared_ptr<SpatialParams> spatialParams_;

    std::shared_ptr<VtkMultiWriter> resultWriter_;
};
} // namespace Dumux

#include <dumux/implicit/adaptive/gridadaptpropertydefaults.hh>

#endif
