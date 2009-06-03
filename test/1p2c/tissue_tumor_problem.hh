// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#ifndef DUNE_TISSUE_TUMOR_PROBLEM_HH
#define DUNE_TISSUE_TUMOR_PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

//#include <dune/grid/alugrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include <dumux/material/phaseproperties/phaseproperties1p.hh>
#include <dumux/boxmodels/1p2c/1p2cboxmodel.hh>
#include "tissue_soilproperties.hh"

/**
 * @file
 * @brief  Definition of a problem, where the distribution of a therapeutic agent
 * within pulmonary tissue is described
 * @author Karin Erbertseder, Bernd Flemisch
 */
namespace Dune
{
template <class TypeTag>
class TissueTumorProblem;

namespace Properties
{
NEW_TYPE_TAG(TissueTumorProblem, INHERITS_FROM(BoxOnePTwoC));

SET_PROP(TissueTumorProblem, Grid)
{
    typedef Dune::UGGrid<2> type;
    //typedef Dune::YaspGrid<2> type;
};

SET_TYPE_PROP(TissueTumorProblem, Problem, Dune::TissueTumorProblem<TTAG(TissueTumorProblem)>);
}

template <class TypeTag = TTAG(TissueTumorProblem) >
class TissueTumorProblem : public BasicDomain<typename GET_PROP_TYPE(TypeTag, PTAG(Grid)),
                                                 typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))      Model;
    typedef typename GridView::Grid                           Grid;

    typedef BasicDomain<Grid, Scalar>    ParentType;
    typedef TissueTumorProblem<TypeTag>      ThisType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // indices of the primary variables
        konti      		 = Indices::konti,
        transport        = Indices::transport,

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

    typedef Dune::InterstitialFluid        Fluid;
    typedef Dune::TissueSoil<Grid, Scalar> Soil;

public:
    TissueTumorProblem(Grid *grid,
                        Scalar dtInitial,
                        Scalar tEnd)
        : ParentType(grid),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("new1p2c")
    {
        timeManager_.setStepSize(dtInitial);

        // specify the grid dimensions
        //outerLowerLeft_[0] = 0;
        //outerLowerLeft_[1] = 0;
        //outerUpperRight_[0] = 50.0;
        //outerUpperRight_[1] = 40.0;

        //height_ = outerUpperRight_[1] - outerLowerLeft_[1];
        //width_  = outerUpperRight_[0] - outerLowerLeft_[0];
        eps_    = 1e-8;

        // for defining e.g. a lense
        //innerLowerLeft_[0] = 0.0;
        //innerLowerLeft_[1] = 0.0;

        //innerUpperRight_[0] = 0.0;
        //innerUpperRight_[1] = 0.5;

        gravity_ = 0;
        //gravity_[dim - 1] = -9.81;
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
    // Strings pulled by the OnePTwoCBoxModel during the course of
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

    //! properties of the liquid phase
    /*! properties of the liquid phase
    */
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        if (globalPos[0] > 22 - eps_)
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        values = 0;

        //Scalar lambda = (globalPos[1])/height_;
        if (globalPos[0] < eps_ ) {
        	values[konti] = -3.8676e-8;
            values[transport] = -4.35064e-10;
        }
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
    	 const GlobalPosition &globalPos
    	            = element.geometry().corner(scvIdx);

        values = Scalar(0.0);

        if(globalPos[0]>10 && globalPos[0]<12 && globalPos[1]>10 && globalPos[1]<12)
        	values[0]= 1.5e-6;
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    Scalar temperature() const
    {
        return 273.15 + 36.0; // 36°C
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
    // the internal method for the initial condition
    void initial_(PrimaryVarVector       &values,
                  const GlobalPosition   &globalPos) const
    {

        values[konti] = -1067;  			//initial condition for the pressure
        values[transport] = 0;   	//initial condition for the molefraction

     }

    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }

    GlobalPosition outerLowerLeft_;
    GlobalPosition outerUpperRight_;
    GlobalPosition innerLowerLeft_;
    GlobalPosition innerUpperRight_;
    Scalar width_;
    Scalar height_;
    Scalar eps_;
    
    GlobalPosition  gravity_;

    // fluids and material properties
    Fluid		    fluid_;
    Soil            soil_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

}; //end namespace
}
#endif
