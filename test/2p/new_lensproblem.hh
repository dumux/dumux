#ifndef DUNE_NEW_LENSPROBLEM_HH
#define DUNE_NEW_LENSPROBLEM_HH

//#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/phaseproperties/phaseproperties2p.hh>


#include <dumux/auxiliary/timemanager.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/new_models/2p/2pboxmodel.hh>

#include <dumux/nonlinear/new_newtonmethod.hh>
#include <dumux/nonlinear/new_newtoncontroller.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include "lenssoil.hh"

namespace Dune
{
/*!
 * \todo Please doc me!
 */
template<class GridT, class ScalarT>
class NewLensProblem : public BasicDomain<GridT,
                                          ScalarT>
{
    typedef GridT                          Grid;
    typedef BasicDomain<Grid, ScalarT>     ParentType;
    typedef NewLensProblem<Grid, ScalarT>  ThisType;
    typedef TwoPBoxModel<ThisType>         Model;

    typedef Water                    WettingPhase;
    typedef DNAPL                    NonwettingPhase;
    typedef Dune::LensSoil<Grid, ScalarT>  Soil;
    typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the Pw-Sn model
    typedef typename Model::TwoPTraits          TwoPTraits;

private:
    // some constants from the traits for convenience
    enum {
        numEq       = BoxTraits::numEq,
        pressureIdx = TwoPTraits::pressureIdx,
        switchIdx   = TwoPTraits::saturationIdx,

        pWIdx       = TwoPTraits::pressureIdx,
        sNIdx       = TwoPTraits::saturationIdx,

        // Grid and world dimension
        dim         = DomainTraits::dim,
        dimWorld    = DomainTraits::dimWorld,
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
    typedef Dune::NewtonController<NewtonMethod>        NewtonController;

public:
    NewLensProblem(Grid *grid,
                   const GlobalPosition &outerLowerLeft,
                   const GlobalPosition &outerUpperRight,
                   const GlobalPosition &innerLowerLeft,
                   const GlobalPosition &innerUpperRight,
                   Scalar dtInitial,
                   Scalar tEnd)
        : ParentType(grid),

          outerLowerLeft_(outerLowerLeft),
          outerUpperRight_(outerUpperRight),

          soil_(outerLowerLeft, outerUpperRight, innerLowerLeft, innerUpperRight),

          materialLaw_(soil_, wPhase_, nPhase_),
          timeManager_(tEnd,
                       this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("newlens")
    {
        timeManager_.setStepSize(dtInitial);

        eps_    = 1e-8 * 300.0;

        gravity_ = 0;
        gravity_[dim - 1] = -9.81;
        
        wasRestarted_ = false;
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

        if (!wasRestarted_) {
            // write the inital solution to disk
            writeCurrentResult_();
        }
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
        model_.update(stepSize,
                      nextStepSize,
                      newtonMethod_,
                      newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        if (this->grid().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        // write the current result to disk
        writeCurrentResult_();

        // write restart file after every five steps
        static int dummy = 0;
        ++dummy;
        if (dummy % 5 == 0)
            serialize();
    };
    ///////////////////////////////////
    // End of simulation control stuff
    ///////////////////////////////////

    ///////////////////////////////////
    // Strings pulled by the TwoPBoxModel during the course of
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

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    MaterialLaw &materialLaw ()
    { return materialLaw_; }

    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = BoundaryConditions::neumann;

        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            values = BoundaryConditions::dirichlet;
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
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);
        
        Scalar densityW = wettingPhase().density();
        
        if (onLeftBoundary_(globalPos))
        {
            Scalar height = outerUpperRight_[1] - outerLowerLeft_[1];
            
            Scalar a = -(1 + 0.5/height);
            Scalar b = -a*outerUpperRight_[1];
            values[pWIdx] = -densityW*gravity_[1]*(a*globalPos[1] + b);
            values[sNIdx] = 0;
        }
        else if (onRightBoundary_(globalPos))
        {
            Scalar a = -1;
            Scalar b = outerUpperRight_[1];
            values[pWIdx] = -densityW*gravity_[1]*(a*globalPos[1] + b);
            values[sNIdx] = 0;
        }
        else
            values = 0.0;
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
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = 0.0;
        /* if (onInlet_(globalPos)) {
            values[sNIdx] = -0.04; // kg / (m * s)
        }
        */
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(SolutionVector &values,
                const Element &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        values = Scalar(0.0);
    }

    //////////////////////////////

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    void initial(SolutionVector         &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);


        Scalar densityW = wettingPhase().density();
        Scalar height = outerUpperRight_[1] - outerLowerLeft_[1];
        values[pWIdx] = -densityW*gravity_[1]*(height - globalPos[1]);

        if (!onLeftBoundary_(globalPos)) {
            Scalar a = -(1 + 0.5/height);
            Scalar b = -a*outerUpperRight_[1];
            values[pWIdx] = -densityW*gravity_[1]*(a*globalPos[1] + b);
        }
        
        values[sNIdx] = 0.0;
    }

    Scalar temperature() const
    {
        return 283.15; // -> 10°C
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

    Model &model()
    {
        return model_;
    }

    const Model &model() const
    {
        return model_;
    }

    void serialize()
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.serializeBegin(this->grid(),
                           "newlens",
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

        res.serializeEnd();
    }

    void deserialize(double t)
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.deserializeBegin(this->grid(), "newlens", t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < outerLowerLeft_[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > outerUpperRight_[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < outerLowerLeft_[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > outerUpperRight_[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = outerUpperRight_[0] - outerLowerLeft_[0];
        Scalar lambda = (outerUpperRight_[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda  && lambda < 2.0/3.0;
    }
    
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }


    Scalar eps_;
    GlobalPosition  gravity_;

    GlobalPosition outerLowerLeft_;
    GlobalPosition outerUpperRight_;

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

    bool wasRestarted_;
};
} //end namespace

#endif
