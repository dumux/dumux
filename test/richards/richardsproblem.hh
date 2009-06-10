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
#ifndef DUNE_RICHARDSPROBLEM_HH
#define DUNE_RICHARDSPROBLEM_HH

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

#include <dune/common/timer.hh>

#include <dumux/material/phaseproperties/phaseproperties2p.hh>

#include<dumux/boxmodels/richards/richardsboxmodel.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include "richardssoil.hh"

namespace Dune
{
/*!
 * \ingroup RichardsBoxProblems
 * \brief  Base class for defining an instance of a Richard`s problem, where water is inflitrating into an initially unsaturated zone.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary) except the top boundary(Dirichlet),
 * where water is inflitrating into an initially unsaturated zone. Linear capillary pressure-saturation relationship is used.
 *
 * To run the simulation execute the following line in shell:
 * ./new_test_richards ./grids/richards_2d.dgf 400000 1
 * where start simulation time = 1 second, end simulation time = 400000 seconds
 * The same file can be also used for 3d simulation but you need to change line
 * "typedef Dune::UGGrid<2> type;" with "typedef Dune::UGGrid<3> type;" and use richards_3d.dgf grid
 */
template <class TypeTag>
class RichardsProblem;

namespace Properties
{
NEW_TYPE_TAG(RichardsProblem, INHERITS_FROM(BoxRichards));

SET_PROP(RichardsProblem, Grid)
{
    typedef Dune::UGGrid<2> type;
//    typedef Dune::ALUCubeGrid<3,3> type;
};

SET_TYPE_PROP(RichardsProblem, Problem, Dune::RichardsProblem<TTAG(RichardsProblem)>);
}

template <class TypeTag = TTAG(RichardsProblem) >
class RichardsProblem : public BasicDomain<typename GET_PROP_TYPE(TypeTag, PTAG(Grid)),
                                              typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))      Model;
    typedef typename GridView::Grid                           Grid;

    typedef BasicDomain<Grid, Scalar>    ParentType;
    typedef RichardsProblem<TypeTag>      ThisType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        pWIdx       = Indices::pWIdx,

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

    typedef Dune::Water                            WettingPhase;
    typedef Dune::DNAPL                            NonwettingPhase;
    typedef Dune::RichardsSoil<Grid, Scalar>       Soil;
    typedef Dune::TwoPhaseRelations<Grid, Scalar>  MaterialLaw;

public:
    RichardsProblem(Grid *grid,
                       Scalar dtInitial,
                       Scalar tEnd)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("richards")
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
    // Strings pulled by the TwoPTwoCBoxModel during the course of
    // the simulation (-> boundary conditions, initial conditions,
    // etc)
    ///////////////////////////////////
    //! Returns the current time step size in seconds
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    //! Set the time step size in seconds.
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }


    //! properties of the wetting (liquid) phase
    /*! properties of the wetting (liquid) phase
      \return    wetting phase
    */
    const WettingPhase &wettingPhase() const
    { return wPhase_; }

    //! properties of the nonwetting (liquid) phase
    /*! properties of the nonwetting (liquid) phase
      \return    nonwetting phase
    */
    const NonwettingPhase &nonwettingPhase() const
    { return nPhase_; }


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

    Model &model()
    {
        return model_;
    }

    const Model &model() const
    {
        return model_;
    }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    MaterialLaw &materialLaw ()
    //        const MaterialLaw &materialLaw () const
    {
        return materialLaw_;
    }

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

        if (globalPos[dim-1] > 20 - eps)
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

         if (globalPos[dim-1] > 20 -eps)
         {
        	values[pWIdx] = 1.0e+5*0.99;
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
//      const GlobalPosition &globalPos = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
    	values[pWIdx] = 0;
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
//        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values[pWIdx] = 1.0e+5 - 5.0e+4;
    }


    Scalar temperature() const
    {
        return 283.15; // 10°C
    };

    Scalar pNreference() const
    {
        return 1.0e+5; // reference non-wetting phase pressure [Pa] used for viscosity and density calculations
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
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
};
} //end namespace

#endif
