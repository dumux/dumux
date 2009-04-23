// $Id$

#ifndef DUNE_NEW_TISSUE_TUMOR_PROBLEM_HH
#define DUNE_NEW_TISSUE_TUMOR_PROBLEM_HH



#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include <dumux/material/phaseproperties/phaseproperties1p.hh>
#include "tissue_soilproperties.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/istl/io.hh>

#include<dumux/new_models/1p2c/1p2cboxmodel.hh>
#include<dumux/new_models/1p2c/1p2cnewtoncontroller.hh>

#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>


/**
 * @file
 * @brief  Definition of a problem, where the distribution of a therapeutic agent
 * within pulmonary tissue is described
 * @author Bernd Flemisch, Karin erbertseder
 */
namespace Dune
{

template<class GridT, class ScalarT>
class NewTissueTumorProblem : public BasicDomain<GridT,
                                               ScalarT>
{
    typedef GridT                               Grid;
    typedef BasicDomain<Grid, ScalarT>          ParentType;
    typedef NewTissueTumorProblem<GridT, ScalarT> ThisType;
    typedef OnePTwoCBoxModel<ThisType>          Model;

    typedef Dune::InterstitialFluid             Phase;
    typedef Dune::TissueSoil<Grid, ScalarT>     Soil;


public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the model
    typedef typename Model::OnePTwoCTraits      OnePTwoCTraits;

private:
    // some constants from the traits for convenience
    enum {
        numEq       = BoxTraits::numEq,

        // Grid and world dimension
        dim          = DomainTraits::dim,
        dimWorld     = DomainTraits::dimWorld,

        // Choice of primary variables
        konti      		 = OnePTwoCTraits::konti,
        transport        = OnePTwoCTraits::transport
    };

    // copy some types from the traits for convenience
    typedef typename DomainTraits::Scalar                     Scalar;
    typedef typename DomainTraits::Element                    Element;
    typedef typename DomainTraits::ElementIterator            ElementIterator;
    typedef typename DomainTraits::ReferenceElement           ReferenceElement;
    typedef typename DomainTraits::Vertex                     Vertex;
    typedef typename DomainTraits::VertexIterator             VertexIterator;
    typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
    typedef typename DomainTraits::LocalPosition              LocalPosition;
    typedef typename DomainTraits::GlobalPosition             GlobalPosition;

    typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction               SpatialFunction;
    typedef typename BoxTraits::SolutionVector                SolutionVector;
    typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

    typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>                  TimeManager;

    typedef typename Model::NewtonMethod                NewtonMethod;
    typedef OnePTwoCNewtonController<NewtonMethod>      NewtonController;

public:
    NewTissueTumorProblem(Grid *grid,
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
        //                writeCurrentResult_(); // TODO
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

        // update the domain with the current solution
        //                updateDomain_();
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




    //! properties of the liquid phase
    /*! properties of the liquid phase
    */
    const Phase &phase() const
    { return phase_; }


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
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        if (globalPos[0] > 22 - eps_)
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    void dirichlet(SolutionVector             &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        initial_(values, globalPos);
    }

    /////////////////////////////
    // NEUMANN boundaries
    /////////////////////////////
    void neumann(SolutionVector             &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);
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
    void source(SolutionVector          &values,
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
    void initial(SolutionVector          &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        /*                const LocalPosition &localPos
                          = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);
        */
        initial_(values, globalPos);
    }


    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };


private:
    // the internal method for the initial condition
    void initial_(SolutionVector         &values,
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


    // fluids and material properties
    Phase		    phase_;
    Soil            soil_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

}; //end namespace
}
#endif
