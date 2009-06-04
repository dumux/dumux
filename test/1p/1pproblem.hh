/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
#ifndef DUNE_1PPROBLEM_HH
#define DUNE_1PPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include <dune/grid/alugrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>

#include<dumux/boxmodels/1p/1pboxmodel.hh>

#include<dumux/nonlinear/newtonmethod.hh>
#include<dumux/nonlinear/newtoncontroller.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include "1psoil.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a single phase isothermal problem
 * @author Onur Dogan, Andreas Lauser
 */
namespace Dune
{
//! base class that defines the parameters of the single phase flow equation
/*!
 * The domain is box shaped. All sides are closed (Neumann 0 boundary) except the top and bottom boundaries,
 * where air is flowing from bottom to top.
 *
 * To run the simulation execute the following line in shell:
 * ./test_1p ./grids/1p_2d.dgf 4000000 1
 * where start simulation time = 1 second, end simulation time = 4000000 seconds
 * The same file can be also used for 3d simulation but you need to change line
 * "typedef Dune::UGGrid<2> type;" with "typedef Dune::UGGrid<3> type;" and use 1p_3d.dgf grid
 */
template <class TypeTag>
class OnePProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePProblem, INHERITS_FROM(BoxOneP));

SET_PROP(OnePProblem, Grid)
{
    typedef Dune::UGGrid<2> type;
//    typedef Dune::ALUCubeGrid<3,3> type;
};

SET_TYPE_PROP(OnePProblem, Problem, Dune::OnePProblem<TTAG(OnePProblem)>);
}


/*!
 * \todo Please doc me!
 */
template <class TypeTag = TTAG(OnePProblem) >
class OnePProblem : public BasicDomain<typename GET_PROP_TYPE(TypeTag, PTAG(Grid)),
                                          typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))      Model;
    typedef typename GridView::Grid                           Grid;

    typedef BasicDomain<Grid, Scalar>    ParentType;
    typedef OnePProblem<TypeTag>      ThisType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        pressureIdx = Indices::pressureIdx,

        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;

    typedef typename GridView::template Codim<0>::Entity    Element;
    typedef typename GridView::template Codim<dim>::Entity  Vertex;
    typedef typename GridView::IntersectionIterator         IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>           TimeManager;
    typedef Dune::VtkMultiWriter<GridView>       VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))      NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController))  NewtonController;

    typedef Dune::Air                     Fluid;
    typedef Dune::OnePSoil<Grid, Scalar>  Soil;

public:
    OnePProblem(Grid *grid,
                   Scalar dtInitial,
                   Scalar tEnd)
        : ParentType(grid),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("1p")
    {
        timeManager_.setStepSize(dtInitial);

        gravity_ = 0;
        gravity_[dim - 1] = -9.81;
    }

    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // set the initial condition
        model_.initial();

        // write the inital solution to disk
        writeCurrentResult_();
    }

    /*!
     * \brief Called by the TimeManager in order to get a time
     *        integration on the model.
     *
     * \note timeStepSize and nextStepSize are references and may
     *       be modified by the TimeIntegration. On exit of this
     *       function 'timeStepSize' must contain the step size
     *       actually used by the time integration for the current
     *       steo, and 'nextStepSize' must contain the suggested
     *       step size for the next time step.
     */
    void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
    {
        model_.update(stepSize, nextStepSize, newtonMethod_, newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        if (this->grid().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        // write the current result to disk
        writeCurrentResult_();
    };
    ///////////////////////////////////
    // End of simulation control stuff
    ///////////////////////////////////

    ///////////////////////////////////
    // Strings pulled by the OnePBoxModel during the course of
    // the simulation (-> boundary conditions, initial conditions,
    // etc)
    ///////////////////////////////////
    //! Returns the current time step size in seconds
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    //! Set the time step size in seconds.
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }

    Model &model()
    {
        return model_;
    }

    const Model &model() const
    {
        return model_;
    }

    //! properties of the fluid
    const Fluid &fluid() const
    { return fluid_; }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    const Soil &soil() const
    {  return soil_; }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    Soil &soil()
    {  return soil_; }

    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
    	double eps = 1.0e-3;
        const GlobalPosition &globalPos
            = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

        if (globalPos[dim-1] < eps || globalPos[dim-1] > 20 - eps)
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
    	double eps = 1.0e-3;
    	const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        if (globalPos[dim-1] < eps) {
            values[pressureIdx] = 2.0e+5;
        }
        else if (globalPos[dim-1] > 20 -eps) {
        	values[pressureIdx] = 1.0e+5;
        }
    }

    /////////////////////////////
    // NEUMANN boundaries
    /////////////////////////////
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
    	//  const GlobalPosition &globalPos = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

    	values[pressureIdx] = 0;
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
        values = Scalar(0.0);
    }

    //////////////////////////////

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values[pressureIdx] = 1.0e+5 + 9.81*1.23*(20-globalPos[dim-1]);
    }

    Scalar temperature() const
    {
        return 283.15; // 10°C
    };

    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };


private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }

    GlobalPosition  gravity_;

    // fluids and material properties
    Fluid           fluid_;
    Soil            soil_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
};
} //end namespace

#endif
