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
#ifndef DUMUX_IMPLICIT_PROBLEM_HH
#define DUMUX_IMPLICIT_PROBLEM_HH

#include "implicitproperties.hh"
#include "implicitmodel.hh"

#include <dumux/io/restart.hh>
#include <dumux/implicit/adaptive/gridadapt.hh>
#include <dumux/common/boundingboxtree.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitBaseProblems
 * \brief Base class for all fully implicit problems
 *
 * \note All quantities are specified assuming a threedimensional
 *       world. Problems discretized using 2D grids are assumed to be
 *       extruded by \f$1 m\f$ and 1D grids are assumed to have a
 *       cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class ImplicitProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef typename GET_PROP_TYPE(TypeTag, VtkMultiWriter) VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    typedef ImplicitGridAdapt<TypeTag, adaptiveGrid> GridAdaptModel;
    typedef BoundingBoxTree<GridView> BoundingBoxTree;

    // copying a problem is not a good idea
    ImplicitProblem(const ImplicitProblem &);

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    ImplicitProblem(TimeManager &timeManager, const GridView &gridView)
        : gridView_(gridView)
        , bBoxMin_(std::numeric_limits<double>::max())
        , bBoxMax_(-std::numeric_limits<double>::max())
        , elementMapper_(gridView)
        , vertexMapper_(gridView)
        , timeManager_(&timeManager)
        , newtonMethod_(asImp_())
        , newtonCtl_(asImp_())
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : Dune::vertices(gridView)) {
            for (int i=0; i<dimWorld; i++) {
                bBoxMin_[i] = std::min(bBoxMin_[i], vertex.geometry().corner(0)[i]);
                bBoxMax_[i] = std::max(bBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        if (gridView.comm().size() > 1)
            for (int i = 0; i < dimWorld; ++i) {
                bBoxMin_[i] = gridView.comm().min(bBoxMin_[i]);
                bBoxMax_[i] = gridView.comm().max(bBoxMax_[i]);
            }

        // set a default name for the problem
        simName_ = "sim";

        // if we are calculating on an adaptive grid get the grid adapt model
        if (adaptiveGrid)
            gridAdapt_ = std::make_shared<GridAdaptModel>(asImp_());
    }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // set the initial condition of the model
        model().init(asImp_());

        // initialize grid adapt model if we have an adaptive grid
        if (adaptiveGrid)
        {
            gridAdapt().init();
        }

        // get and apply point sources if any given in the problem
        std::vector<PointSource> sources;
        asImp_().addPointSources(sources);
        // if there are point sources compute the DOF to point source map
        if (!sources.empty())
        {
            // build the bounding box tree for fast point in element search
            boundingBoxTree_ = std::make_shared<BoundingBoxTree>(gridView_);

            // calculate point source locations and save them in a map
            Dumux::PointSourceHelper<TypeTag>::computePointSourceMap(asImp_(),
                                                                     boundingBoxTree_,
                                                                     sources,
                                                                     pointSourceMap_);
        }
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values,
                       const Vertex &vertex) const
    {
        if (!isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "boundaryTypes(..., vertex) called for cell-centered method.");

        // forward it to the method which only takes the global coordinate
        asImp_().boundaryTypesAtPos(values, vertex.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param intersection The intersection for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values,
                       const Intersection &intersection) const
    {
        if (isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "boundaryTypes(..., intersection) called for box method.");

        // forward it to the method which only takes the global coordinate
        asImp_().boundaryTypesAtPos(values, intersection.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the finite volume in global coordinates
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a boundaryTypes() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values,
                   const Vertex &vertex) const
    {
        if (!isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "dirichlet(..., vertex) called for cell-centered method.");

        // forward it to the method which only takes the global coordinate
        asImp_().dirichletAtPos(values, vertex.geometry().center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param intersection The intersection for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values,
                   const Intersection &intersection) const
    {
        if (isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "dirichlet(..., intersection) called for box method.");

        // forward it to the method which only takes the global coordinate
        asImp_().dirichletAtPos(values, intersection.geometry().center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
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
                   "a dirichlet() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        // forward it to the interface without the volume variables
        asImp_().neumann(values,
                         element,
                         fvGeometry,
                         intersection,
                         scvIdx,
                         boundaryFaceIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        // forward it to the interface with only the global position
        asImp_().neumannAtPos(values, fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param globalPos The position of the boundary face's integration point in global coordinates
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
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void solDependentSource(PrimaryVariables &values,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int scvIdx,
                     const ElementVolumeVariables &elemVolVars) const
    {
        // forward to solution independent, fully-implicit specific interface
        asImp_().source(values, element, fvGeometry, scvIdx);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx) const
    {
        // forward to generic interface
        asImp_().sourceAtPos(values, fvGeometry.subContVol[scvIdx].global);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
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
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a sourceAtPos() method.");
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const {}

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSources The vector of point sources
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void solDependentPointSources(std::vector<PointSource> &pointSources,
                                  const Element &element,
                                  const FVElementGeometry &fvGeometry,
                                  const int scvIdx,
                                  const ElementVolumeVariables &elemVolVars) const {}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        // forward to generic interface
        asImp_().initialAtPos(values, fvGeometry.subContVol[scvIdx].global);
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
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a initialAtPos() method.");
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
    Scalar boxExtrusionFactor(const Element &element,
                              const FVElementGeometry &fvGeometry,
                              const int scvIdx) const
    {
        // forward to generic interface
        return asImp_().extrusionFactorAtPos(fvGeometry.subContVol[scvIdx].global);
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
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    { return 1.0; }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void preTimeStep()
    {
        // If adaptivity is used, this method adapts the grid.
        // Remeber to call the parent class function if this is overwritten
        // on a lower problem level when using an adaptive grid
        if (adaptiveGrid && timeManager().timeStepIndex() > 0)
            this->gridAdapt().adaptGrid();
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
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

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " time-step divisions. dt="
                   << timeManager().timeStepSize());
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
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        time step has been computed and the simulation time has
     *        been updated.
     *
     * \param dt The current time-step size
     */
    Scalar nextTimeStepSize(const Scalar dt)
    {
        return newtonCtl_.suggestTimeStepSize(dt);
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 5 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        return timeManager().timeStepIndex() > 0 &&
            (timeManager().timeStepIndex() % 10 == 0);
    }

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
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    { }

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
    const char *name() const
    {
        return simName_.c_str();
    }

    /*!
     * \brief Set the problem name.
     *
     * This static method sets the simulation name, which should be
     * called before the application problem is declared! If not, the
     * default name "sim" will be used.
     *
     * \param newName The problem's name
     */
    void setName(const char *newName)
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
     * file.)  See Dumux::Restart for details.
     */
    void serialize()
    {
        typedef Dumux::Restart Restarter;
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
     * file.)  See Dumux::Restart for details.
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
        typedef Dumux::Restart Restarter;

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
        if (asImp_().shouldWriteOutput()) {
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
    }

    /*!
     * \brief Returns a reference to the grid
     */
    Grid &grid()
    {
        return GridCreator::grid();
    }

    /*!
     * \brief Returns adaptivity model used for the problem.
     */
    GridAdaptModel& gridAdapt()
    {
        return *gridAdapt_;
    }

    /*!
     * \brief Returns adaptivity model used for the problem.
     */
    const GridAdaptModel& gridAdapt() const
    {
        return *gridAdapt_;
    }

    /*!
     * \brief Returns the bounding box tree of the grid
     */
    BoundingBoxTree& boundingBoxTree()
    {
        if(!boundingBoxTree_)
            boundingBoxTree_ = std::make_shared<BoundingBoxTree>(gridView_);

        return *boundingBoxTree_;
    }

    /*!
     * \brief Returns the bounding box tree of the grid
     */
    const BoundingBoxTree& boundingBoxTree() const
    {
        if(!boundingBoxTree_)
            boundingBoxTree_ = std::make_shared<BoundingBoxTree>(gridView_);

        return *boundingBoxTree_;
    }

    /*!
     * \brief Adds contribution of point sources for a specific DOF
     *        to the values.
     */
    void dofSources(PrimaryVariables &values,
                    const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx,
                    const ElementVolumeVariables &elemVolVars) const
    {
        unsigned int dofGlobalIdx = model_.dofMapper().subIndex(element, scvIdx, dofCodim);
        if (pointSourceMap_.count(dofGlobalIdx))
        {
            // call the solDependent function. Herein the user might fill/add values to the point sources
            auto pointSources = pointSourceMap_.at(dofGlobalIdx);
            asImp_().solDependentPointSources(pointSources, element, fvGeometry, scvIdx, elemVolVars);

            // Add the contributions to the dof source values
            Scalar volume = isBox ? model_.boxVolume(dofGlobalIdx) : element.geometry().volume();
            for (auto&& pointSource : pointSources)
            {
                pointSource /= volume;
                values += pointSource.values();
            }
        }
    }

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
    //! \copydoc Dumux::IMPETProblem::resultWriter()
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

        // Tell the result writer that the grid changes if we are adaptive
        if (adaptiveGrid)
        {
            resultWriter_->gridChanged();
        }
    }

    std::string simName_;
    const GridView gridView_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    TimeManager *timeManager_;

    Model model_;

    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;

    std::shared_ptr<VtkMultiWriter> resultWriter_;

    std::shared_ptr<GridAdaptModel> gridAdapt_;

    std::shared_ptr<BoundingBoxTree> boundingBoxTree_;
    std::map<unsigned int, std::vector<PointSource> > pointSourceMap_;
};
} // namespace Dumux

#include <dumux/implicit/adaptive/gridadaptpropertydefaults.hh>

#endif
