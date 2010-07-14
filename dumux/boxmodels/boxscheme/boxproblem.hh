// $Id: boxproblem.hh 3784 2010-06-24 13:43:57Z bernd $
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

#include <dumux/common/timemanager.hh>

namespace Dumux
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
    typedef Dumux::TimeManager<Episode>      TimeManager;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))      NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController))  NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))             Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))            Scalar;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::Grid::ctype                         CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld>             GlobalPosition;
    typedef typename GridView::template Codim<dim>::Iterator     VertexIterator;

public:
    BoxProblem(const GridView &gridView)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          resultWriter_(asImp_()->name())
    {
        wasRestarted_ = false;

        // calculate the bounding box of the grid view
        VertexIterator vIt = gridView.template begin<dim>();
        const VertexIterator vEndIt = gridView.template end<dim>();
        for (; vIt!=vEndIt; ++vIt) {
            for (int i=0; i<dim; i++) {
                bboxMin_[i] = std::min(bboxMin_[i], vIt->geometry().corner(0)[i]);
                bboxMax_[i] = std::max(bboxMax_[i], vIt->geometry().corner(0)[i]);
            }
        }
        for (int i=0; i<dim; i++) {
            bboxMin_[i] = gridView.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView.comm().max(bboxMax_[i]);
        }
        model_ = new Model(*asImp_());
        newtonMethod_ = new NewtonMethod(*model_);
    }

    ~BoxProblem()
    {
        delete model_;
        delete newtonMethod_;
    };

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Start the simulation procedure.
     *
     * This method is usually called by the main() function and simply
     * uses Dumux::TimeManager::runSimulation() to do the actual
     * work.
     */
    bool simulate(Scalar dtInitial, Scalar tEnd)
    {
        // set the initial time step and the time where the simulation ends
        timeManager_.setEndTime(tEnd);
        timeManager_.setTimeStepSize(dtInitial);
        timeManager_.runSimulation(*asImp_());
        return true;
    };


    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model().initial();

        // write the inital solution to disk
        asImp_()->writeCurrentResult_();
    }

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void timeStepBegin()
    {}

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        model().update(*newtonMethod_, newtonCtl_);
    }

    /*!
     * \brief Called by Dumux::TimeManager whenever a solution for a
     *        timestep has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize()
    {
        Scalar dt = asImp_()->timeManager().timeStepSize();
        return newtonCtl_.suggestTimeStepSize(dt);
    };


    /*!
     * \brief This method is called by the model if the update to the
     *        next time step failed completely.
     */
    void updateSuccessful()
    {
        wasRestarted_ = false;
        asImp_()->writeCurrentResult_();
    };

    /*!
     * \brief This method is called by the model if the update to the
     *        next time step failed completely.
     */
    void updateFailed()
    { };

    /*!
     * \brief This method is called by the model if the update to the
     *        next time step failed with the curent time step size.
     */
    void updateFailedTry()
    {
        Scalar dt = asImp_()->timeManager().timeStepSize();
        Scalar nextDt = newtonCtl_.suggestTimeStepSize(dt);
        asImp_()->timeManager().setTimeStepSize(nextDt);
    };

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
            timeManager().timeStepNum() > 0 &&
            (timeManager().timeStepNum() % 5 == 0);
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

    /*!
     * \brief Called by the time manager after the time integration.
     */
    void timeStepEnd()
    {
        // write restart file if necessary
        if (asImp_()->shouldWriteRestartFile())
            serialize();
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

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model &model()
    { return *model_; }

    /*!
     * \copydoc model()
     */
    const Model &model() const
    { return *model_; }
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
     * file.)  See Dumux::Restart for details.
     */
    void serialize()
    {
        typedef Dumux::Restart<GridView> Restarter;

        Restarter res;
        res.serializeBegin(gridView(),
                           asImp_()->name(),
                           timeManager_.time());
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
        typedef Dumux::Restart<GridView> Restarter;

        Restarter res;
        res.deserializeBegin(gridView(),
                             asImp_()->name(),
                             t);
        std::cerr << "Deserialize from file " << res.fileName() << "\n";

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model().deserialize(res);

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

            // calculate the time _after_ the time was updated
            Scalar t = timeManager_.time() + timeManager_.timeStepSize();
            resultWriter_.beginTimestep(t,
                                        gridView());
            model().addOutputVtkFields(resultWriter_);
            resultWriter_.endTimestep();
        }
    }

private:
    static std::string simname_; // a string for the name of the current simulation,
                                  // which could be set by means of an program argument,
                                  // for example.
    const GridView  gridView_;

    GlobalPosition  bboxMin_;
    GlobalPosition  bboxMax_;

    TimeManager     timeManager_;

    Model           *model_;

    NewtonMethod    *newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

    bool wasRestarted_;
};
// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag, class Implementation>
std::string BoxProblem<TypeTag, Implementation>::simname_="sim"; //initialized with default "sim"

}

#endif
