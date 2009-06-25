// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
/*!
 * \file
 * \brief Base class for all problems which use the box scheme
 */
#ifndef DUMUX_BOX_PROBLEM_HH
#define DUMUX_BOX_PROBLEM_HH

#include "boxproperties.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/auxiliary/timemanager.hh>

namespace Dune
{
/*!
 * \ingroup BoxScheme
 * \brief  Base class for all problems which use the box scheme
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation>
class BoxProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>      TimeManager;

    typedef Dune::VtkMultiWriter<GridView>  VtkMultiWriter;
    
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))      NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController))  NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))             Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))            Scalar;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

public:
    BoxProblem(const GridView &gridView)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          timeManager_(gridView.comm().rank() == 0),
          model_(*asImp_()),
          newtonMethod_(model_),
          resultWriter_(asImp_()->name())
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
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Start the simulation procedure. 
     *
     * This method is usually called by the main() function and simply
     * uses Dune::TimeManager::runSimulation() to do the actual
     * work.
     */
    bool simulate(Scalar dtInitial, Scalar tEnd)
    {
        // set the initial time step and the time where the simulation ends
        timeManager_.setEndTime(tEnd);
        timeManager_.setStepSize(dtInitial);
        timeManager_.runSimulation(*asImp_());
        return true;
    };


    /*!
     * \brief Called by the Dune::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model_.initial();

        // write the inital solution to disk
        writeCurrentResult_();
    }

    /*!
     * \brief Called by Dune::TimeManager in order to do a time
     *        integration on the model.
     *
     * \note \a timeStepSize and \a nextStepSize are references and may
     *       be modified by the timeIntegration(). On exit of this
     *       function \a timeStepSize must contain the step size
     *       actually used by the time integration for the current
     *       steo, and \a nextStepSize must contain a suggestion for the 
     *       next time step size.
     */
    void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
    {
        model_.update(stepSize,
                      nextStepSize,
                      newtonMethod_,
                      newtonCtl_);
    }

    /*!
     * \brief Called by Dune::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     *
     * This is used to do some janitorial tasks like writing the
     * current solution to disk.
     */
    void timestepDone()
    {
        asImp_()->writeCurrentResult_();
        wasRestarted_ = false;
    };

    /*!
     * \brief Returns the current time step size [seconds].
     */
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    /*!
     * \brief Sets the current time step size [seconds].
     */
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }

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
        return !restarted() && 
            timeManager().stepNum() > 0 && 
            (timeManager().stepNum() % 5 == 0);
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behaviour is to write out every the solution for
     * very time step. This file is intented to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutputFile() const 
    { return !restarted(); }

    // \}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation
     * and should be overwritten by the problem if these files should
     * not start with "sim".
     */
    const char *name() const
    { return "sim"; }

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
     * \brief Returns true, if the current state of the problem was
     *        loaded from a restart file.
     */
    bool restarted() const
    { return wasRestarted_; }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.drs</tt>. (Dumux ReStart
     * file.)  See Dune::Restart for details.
     */
    void serialize()
    {
        typedef Dune::Restart<GridView> Restarter;

        Restarter res;
        res.serializeBegin(gridView(),
                           asImp_()->name(),
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

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
        typedef Dune::Restart<GridView> Restarter;

        Restarter res;
        res.deserializeBegin(gridView(),
                             asImp_()->name(),
                             t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };

    // \}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    
    //! \copydoc asImp_()
    const Implementation *asImp_() const 
    { return static_cast<const Implementation *>(this); }

    //! Write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        // write the current result to disk
        if (asImp_()->shouldWriteOutputFile()) {
            if (gridView().comm().rank() == 0)
                std::cout << "Writing result file for current time step\n";

            resultWriter_.beginTimestep(timeManager_.time(),
                                        gridView());
            model_.addVtkFields(resultWriter_);
            resultWriter_.endTimestep();
        }

        // write restart file if necessary
        if (asImp_()->shouldWriteRestartFile())
            serialize();
    }

private:
    const GridView  gridView_;
    
    GlobalPosition  bboxMin_;
    GlobalPosition  bboxMax_;

    TimeManager     timeManager_;

    Model            model_;

    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

    bool wasRestarted_;
};
}

#endif
