// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Andreas Lauser                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_IMPESPROBLEM_HH
#define DUMUX_IMPESPROBLEM_HH

#include "impesproperties.hh"
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/common/timemanager.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the diffusion problem
 * @author Bernd Flemisch
 */

namespace Dumux
{

/*! \ingroup fracflow
 * @brief base class that defines the parameters of loosely coupled diffusion and transport equations
 *
 * An interface for defining parameters for the stationary diffusion equation
 *  \f[\text{div}\, \boldsymbol{v} = q\f]
 *  and a scalar transport equation
 *  \f[
 *    \frac{\partial S}{\partial t} + \text{div}\, \boldsymbol{v_\alpha} = 0,
 *  \f]
 *  where, the velocity \f$\boldsymbol{v} \sim \boldsymbol{K} \nabla p \f$,
 *  \f$p\f$ is a pressure and q a source/sink term, \f$S\f$ denotes a phase saturation and \f$\boldsymbol{v_\alpha}\f$ is a phase velocity.
 *
 *  Template parameters are:
 *
 *  - GridView      a DUNE gridview type
 *  - Scalar        type used for scalar quantities
 *  - VC            type of a class containing different variables of the model
 */
template<class TypeTag, class Implementation>
class IMPESProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef Dumux::TimeManager<TypeTag>  TimeManager;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) IMPESModel;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportSolutionType)) TransportSolutionType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureModel)) PressureModel;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SaturationModel)) SaturationModel;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

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

public:

    //! Constructs an object of type IMPESProblemProblem
    /** @param variables object of class VariableClass.
     *  @param wettingPhase implementation of a wetting phase.
     *  @param nonWettingPhase implementation of a non-wetting phase.
     *  @param spatial parameters implementation of the solid matrix
     *  @param materialLaw implementation of Material laws. Class TwoPhaseRelations or derived.
     */
    IMPESProblem(const GridView &gridView, bool verbose = true)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          timeManager_(verbose),
          variables_(gridView),
          dt_(0),
          resultWriter_(asImp_().name())
    {
        // calculate the bounding box of the grid view
        VertexIterator vIt = gridView.template begin<dim>();
        const VertexIterator vEndIt = gridView.template end<dim>();
        for (; vIt!=vEndIt; ++vIt) {
            for (int i=0; i<dim; i++) {
                bboxMin_[i] = std::min(bboxMin_[i], vIt->geometry().corner(0)[i]);
                bboxMax_[i] = std::max(bboxMax_[i], vIt->geometry().corner(0)[i]);
            }
        }

        pressModel_ = new PressureModel(asImp_());
        satModel_ = new SaturationModel(asImp_());
        model_ = new IMPESModel(asImp_()) ;
    }

    //! destructor
    virtual ~IMPESProblem ()
    {
        delete pressModel_;
        delete satModel_;
        delete model_;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model().initial();
    }

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void preTimeStep()
    {}

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

        // obtain the first update and the time step size
        model().update(t, dt_, k1);

        //make sure t_old + dt is not larger than tend
        dt_ = std::min(dt_, timeManager().episodeMaxTimeStepSize());

        // check if we are in first TS and an initialDt was assigned
        if (t==0. && timeManager().timeStepSize()!=0.)
        {
            // check if assigned initialDt is in accordance with dt_ from first transport step
            if (timeManager().timeStepSize() > dt_)
                Dune::dwarn << "initial timestep of size " << timeManager().timeStepSize()
                            << "is larger then dt_= "<<dt_<<" from transport" << std::endl;
            // internally assign next tiestep size
            dt_ = std::min(dt_, timeManager().timeStepSize());
        }

        //assign next tiestep size
        timeManager().setTimeStepSize(dt_);

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
    Scalar nextTimeStepSize()
    { return dt_;}

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
     * \copydoc timeManager()
     */
    const TimeManager &timeManager() const
    { return timeManager_; }

    Variables& variables ()
    { return variables_; }

    const Variables& variables () const
    { return variables_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    IMPESModel &model()
    { return *model_; }

    /*!
     * \copydoc model()
     */
    const IMPESModel &model() const
    { return *model_; }
    // \}

    /*!
     * \brief Returns numerical model used for the problem.
     */
    PressureModel &pressureModel()
    { return *pressModel_; }

    /*!
     * \copydoc model()
     */
    const PressureModel &pressureModel() const
    { return *pressModel_; }
    // \}

    /*!
     * \brief Returns numerical model used for the problem.
     */
    SaturationModel &saturationModel()
    { return *satModel_; }

    /*!
     * \copydoc model()
     */
    const SaturationModel &saturationModel() const
    { return *satModel_; }
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

    //! Write the fields current solution into an VTK output file.
    void writeOutput()
    {
        if (gridView().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        resultWriter_.beginTimestep(timeManager_.time() + timeManager_.timeStepSize(),
                                    gridView());
        model().addOutputVtkFields(resultWriter_);
        resultWriter_.endTimestep();
    }

    // \}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
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

    Scalar dt_;

    PressureModel* pressModel_;//!< object including the pressure model
    SaturationModel* satModel_;//!< object including the saturation model
    IMPESModel* model_;

    VtkMultiWriter resultWriter_;
};
// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag, class Implementation>
std::string IMPESProblem<TypeTag, Implementation>::simname_="sim"; //initialized with default "sim"
}
#endif
