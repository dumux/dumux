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

#include <dumux/common/timemanager.hh>

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Base class for all problems which use the box scheme
 *
 * \todo Please doc me more!
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

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    // copying a problem is not a good idea
    BoxProblem(const BoxProblem &);

public:
    BoxProblem(TimeManager &timeManager, const GridView &gridView)
        : gridView_(gridView),
          bboxMin_(std::numeric_limits<double>::max()),
          bboxMax_(-std::numeric_limits<double>::max()),
          elementMapper_(gridView),
          vertexMapper_(gridView),
          timeManager_(&timeManager),
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
        for (int i=0; i<dim; i++) {
            bboxMin_[i] = gridView.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView.comm().max(bboxMax_[i]);
        }

        model_ = new Model();
        newtonMethod_ = new NewtonMethod(asImp_());
    }

    ~BoxProblem()
    {
        delete model_;
        delete newtonMethod_;
    };

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model().init(asImp_());
    }

    /*!
     * \brief If model coupling is used, this updates the parameters
     *        required to calculate the coupling fluxes between the
     *        sub-models.
     *
     * By default it does nothing
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

            if (model_->update(*newtonMethod_, newtonCtl_))
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
     */
    Scalar nextTimeStepSize()
    {
        Scalar dt = asImp_().timeManager().timeStepSize();
        return newtonCtl_.suggestTimeStepSize(dt);
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
     * \brief Called by the time manager after the time integration.
     */
    void postTimeStep()
    { }

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
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        resultWriter_.serialize(res);
        model().serialize(res);
    }

    /*!
     * \brief Load a previously saved state of the whole simulation
     *        from disk.
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
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        resultWriter_.deserialize(res);
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
            resultWriter_.beginTimestep(t, gridView());
            model().addOutputVtkFields(model().curSol(), resultWriter_);
            resultWriter_.endTimestep();
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
    static std::string simname_; // a string for the name of the current simulation,
                                  // which could be set by means of an program argument,
                                  // for example.
    const GridView gridView_;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    TimeManager *timeManager_;

    Model *model_;

    NewtonMethod    *newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter resultWriter_;
};
// definition of the static class member simname_,
// which is necessary because it is of type string.
template <class TypeTag>
std::string BoxProblem<TypeTag>::simname_="sim"; //initialized with default "sim"

}

#endif
