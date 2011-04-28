// $Id$
/*****************************************************************************
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
/*!
 * \file
 * \brief Base class for all problems which use the box scheme
 */
#ifndef DUMUX_BOX_PROBLEM_HH
#define DUMUX_BOX_PROBLEM_HH

#include "boxproperties.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Base class for all problems which use the box scheme.
 *
 * \note All quantities are specified assuming a threedimensional
 *       world. Problems discretized using 2D grids are assumed to be
 *       extruded by \f$1 m\f$ and 1D grids are assumed to have a
 *       cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class BoxProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef Dumux::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::Intersection Intersection;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    // copying a problem is not a good idea
    BoxProblem(const BoxProblem &);

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    BoxProblem(TimeManager &timeManager, const GridView &gridView)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          elementMapper_(gridView),
          vertexMapper_(gridView),
          timeManager_(&timeManager),
          newtonMethod_(asImp_())
    {
        // calculate the bounding box of the local partition of the grid view
        VertexIterator vIt = gridView.template begin<dim>();
        const VertexIterator vEndIt = gridView.template end<dim>();
        for (; vIt!=vEndIt; ++vIt) {
            for (int i=0; i<dim; i++) {
                bboxMin_[i] = std::min(bboxMin_[i], vIt->geometry().corner(0)[i]);
                bboxMax_[i] = std::max(bboxMax_[i], vIt->geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        for (int i = 0; i < dim; ++i) {
            bboxMin_[i] = gridView.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView.comm().max(bboxMax_[i]);
        };

        // set a default name for the problem
        simName_ = "sim";

        resultWriter_ = NULL;
    }

    ~BoxProblem()
    {
        delete resultWriter_;
    };


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
    }

    /*!
     * \brief Returns the maximum allowed time step size [s]
     *
     * By default this the time step size is unrestricted.
     */
    Scalar maxTimeStepSize() const 
    { return std::numeric_limits<Scalar>::infinity(); }

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
        // forward it to the method which only takes the global coordinate
        asImp_().boundaryTypes(values, vertex.geometry().center());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param pos The position of the finite volume in global coordinates
     */
    void boundaryTypes(PrimaryVariables &values,
                       const GlobalPosition &pos) const
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
        // forward it to the method which only takes the global coordinate
        asImp_().dirichlet(values, vertex.geometry().center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param pos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values,
                   const GlobalPosition &pos) const
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
     * potentially solution dependent and requires some box method
     * specific things.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void boxSDNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvElemGeom,
                      const Intersection &is,
                      int scvIdx,
                      int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        // forward it to the interface without the volume variables
        asImp_().neumann(values,
                         element,
                         fvElemGeom,
                         is,
                         scvIdx,
                         boundaryFaceIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        // forward it to the interface with only the global position
        asImp_().neumann(values, fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param pos The position of the boundary face's integration point in global coordinates
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const GlobalPosition &pos) const
    {
        // do nothing
        values = 0.0;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some box method
     * specific things.
     *
     * \param values The source and sink values for the conservation equations
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void boxSDSource(PrimaryVariables &values,
                     const Element &element,
                     const FVElementGeometry &fvElemGeom,
                     int scvIdx,
                     const ElementVolumeVariables &elemVolVars) const
    {
        // forward to solution independent, box specific interface
        asImp_().source(values, element, fvElemGeom, scvIdx);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        // forward to generic interface
        asImp_().source(values, fvElemGeom.subContVol[scvIdx].global);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations
     * \param pos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const GlobalPosition &pos) const
    { values = Scalar(0.0);  }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        // forward to generic interface
        asImp_().initial(values, fvElemGeom.subContVol[scvIdx].global);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param pos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void initial(PrimaryVariables &values,
                 const GlobalPosition &pos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a initial() method.");
    }

    /*!
     * \brief If model coupling is used, this updates the parameters
     *        required to calculate the coupling fluxes between the
     *        sub-models.
     *
     * By default it does nothing
     *
     * \param element The DUNE Codim<0> entity for which the coupling
     *                parameters should be computed.
     */
    void updateCouplingParams(const Element &element) const
    {}

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
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        const int maxFails = 10;
        for (int i = 0; i < maxFails; ++i) {
            if (i > 0 && gridView().comm().rank() == 0)
                std::cout << "Newton solver did not converge. Retrying with time step of "
                          << timeManager().timeStepSize() << "sec\n";

            if (model_.update(newtonMethod_, newtonCtl_))
                return;

            // update failed
            Scalar dt = timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            timeManager().setTimeStepSize(nextDt);
        }

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " timestep divisions. dt="
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
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * \param dt The current time step size
     */
    Scalar nextTimeStepSize(Scalar dt)
    { return newtonCtl_.suggestTimeStepSize(dt); };

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
        return timeManager().timeStepIndex() > 0 &&
            (timeManager().timeStepIndex() % 10 == 0);
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
     * \brief Returns the mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
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
        std::cerr << "Serialize to file '" << res.fileName() << "'\n";

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
    void restart(Scalar tRestart)
    {
        typedef Dumux::Restart Restarter;

        Restarter res;

        res.deserializeBegin(asImp_(), tRestart);
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
     * \param res The desrializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        createResultWriter_();
        resultWriter_->deserialize(res);
        model().deserialize(res);
    };

    // \}

    /*!
     * \brief Write the relavant secondar variables of the current
     *        solution into an VTK output file.
     */
    void writeOutput()
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            if (gridView().comm().rank() == 0)
                std::cout << "Writing result file for \"" << asImp_().name() << "\"\n";

            // calculate the time _after_ the time was updated
            Scalar t = timeManager().time() + timeManager().timeStepSize();
            createResultWriter_();
            resultWriter_->beginTimestep(t, gridView());
            model().addOutputVtkFields(model().curSol(), *resultWriter_);
            resultWriter_->endTimestep();
        }
    }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    // makes sure that the result writer exists
    void createResultWriter_()
    { if (!resultWriter_) resultWriter_ = new VtkMultiWriter(asImp_().name()); };

    std::string simName_; 
    const GridView gridView_;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    TimeManager *timeManager_;

    Model model_;

    NewtonMethod newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter *resultWriter_;
};

}

#endif
