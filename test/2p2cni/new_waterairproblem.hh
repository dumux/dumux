#ifndef DUNE_NEW_WATERAIRPROBLEM_HH
#define DUNE_NEW_WATERAIRPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

//#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/multicomponentrelations.hh>

#include <dumux/material/phaseproperties/phaseproperties2p.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/new_models/2p2cni/2p2cniboxmodel.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#define ISOTHERMAL 0

/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
template <class TypeTag>
class NewWaterAirProblem;

namespace Properties
{
NEW_TYPE_TAG(WaterAirProblem, 
#if ISOTHERMAL
             INHERITS_FROM(BoxTwoPTwoC));
#else
             INHERITS_FROM(BoxTwoPTwoCNI));
#endif


SET_PROP(WaterAirProblem, Grid)
{
    typedef Dune::UGGrid<2> type;
};

SET_PROP(WaterAirProblem, Problem)
{
    typedef Dune::NewWaterAirProblem<TTAG(WaterAirProblem)> type;
};
}


//! class that defines the parameters of an air waterair under a low permeable layer
/*! Problem definition of an air injection under a low permeable layer. Air enters the domain
 * at the right boundary and migrates upwards.
 * Problem was set up using the rect2d.dgf grid.
 *
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template <class TypeTag = TTAG(WaterAirProblem) >
class NewWaterAirProblem : public BasicDomain<typename GET_PROP_TYPE(TypeTag, PTAG(Grid)),
                                              typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))      Model;
    typedef typename GridView::Grid                           Grid;

    typedef BasicDomain<Grid, Scalar>    ParentType;
    typedef NewWaterAirProblem<TypeTag>  ThisType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        pressureIdx = Indices::pressureIdx,
        switchIdx   = Indices::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
#endif

        // Phase State
        wPhaseOnly  = Indices::wPhaseOnly,
        nPhaseOnly  = Indices::nPhaseOnly,
        bothPhases  = Indices::bothPhases,

        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef typename GridView::template Codim<dim>::Entity      Vertex;
    typedef typename GridView::IntersectionIterator             IntersectionIterator;
  
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
   

    typedef Dune::Liq_WaterAir                     WettingPhase;
    typedef Dune::Gas_WaterAir                     NonwettingPhase;
    typedef Dune::HomogeneousSoil<Grid, Scalar>    Soil;
    typedef Dune::TwoPhaseRelations<Grid, Scalar>  MaterialLaw;
    typedef Dune::CWaterAir                        Multicomp;

    typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>                  TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))     NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

public:
    NewWaterAirProblem(Grid *grid,
                       Scalar dtInitial,
                       Scalar tEnd)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          multicomp_(wPhase_, nPhase_),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("new_waterair")
    {
        timeManager_.setStepSize(dtInitial);

        // specify the grid dimensions
        width_ = 60; // [m]
        height_ = 40; // [m]

        depthBOR_ = 1000.0;

        gravity_ = 0;
        gravity_[dim-1] = -9.81;
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

    //! object for multicomponent calculations
    /*! object for multicomponent calculations including mass fractions,
     * mole fractions and some basic laws
     \return    multicomponent object
    */
    MultiComp &multicomp ()
    //        const MultiComp &multicomp () const
    {
        return multicomp_;
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        if(globalPos[0] < eps_)
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;

#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
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
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

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
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = 0;

        // negative values for injection
        if (globalPos[0] > width_ - eps_ &&
            globalPos[1] < 5.0 && globalPos[1] < 15.0)
        {
            values[switchIdx] = -1e-3;
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        /*                const LocalPosition &localPos
                          = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);
        */

        initial_(values, globalPos);
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVarVector     &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
        values[switchIdx] = 0.0;
#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
#endif
    }

public:
    int initialPhaseState(const Vertex         &vert,
                          int                  &globalIdx,
                          const GlobalPosition &globalPos) const
    {
        return wPhaseOnly;
    }

#if ISOTHERMAL
    Scalar temperature() const 
    { 
        return 303.0;
    }
#endif

    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    double depthBOR () const
    {
        return depthBOR_;
    }

    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };

    Model &model()
    {
        return model_;
    }

    const Model &model() const
    {
        return model_;
    }

private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }


    Scalar width_;
    Scalar height_;
    Scalar depthBOR_;
    static const Scalar eps_ = 1e-8;
    GlobalPosition  gravity_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;
    Multicomp       multicomp_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
};
} //end namespace

#endif
