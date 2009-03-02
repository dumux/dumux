#ifndef DUNE_TIMEDEPENDENTPROBLEM_HH
#define DUNE_TIMEDEPENDENTPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include <dumux/material/phaseproperties/phaseproperties_waterair.hh>
#include <timeproblem_soilproperties.hh>
#include <dumux/material/twophaserelations.hh>

#include <dumux/material/multicomponentrelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/istl/io.hh>

#include<dumux/new_models/2p2cni/2p2cniboxmodel.hh>
#include<dumux/new_models/2p2cni/2p2cninewtoncontroller.hh>

#include<dumux/timedisc/new_impliciteulerstep.hh>
#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#define ISOTHERMAL 0

/**
 * @file
 * @brief  Definition of a problem with a time dependent boundary condition (artificial)
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
//! class that defines a problem with time dependent boundary conditions
/*! Problem definition, where the lower dirichlet boundary changes over time
 *  The timedepproblem_soil.hh class was used
 *
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template<class ScalarT>
class TimeDependentProblem : public BasicDomain<Dune::SGrid<2,2>,
                                                ScalarT>
{
    typedef Dune::SGrid<2,2>               Grid;
    typedef BasicDomain<Grid, ScalarT>     ParentType;
    typedef TimeDependentProblem<ScalarT>    ThisType;
#if !ISOTHERMAL
    typedef TwoPTwoCNIBoxModel<ThisType>   Model;
#else
    typedef TwoPTwoCBoxModel<ThisType>   Model;
#endif

    typedef Dune::GridPtr<Grid>                    GridPointer;

    typedef Dune::Liq_WaterAir                     WettingPhase;
    typedef Dune::Gas_WaterAir                     NonwettingPhase;
    typedef Dune::TimeDepProblemSoil<Grid, ScalarT>   Soil;
    typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;
    typedef Dune::CWaterAir                        Multicomp;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the Pw-Sn model
#if !ISOTHERMAL
    typedef typename Model::TwoPTwoCNITraits    TwoPTwoCNITraits;
#else
    typedef typename Model::TwoPTwoCTraits      TwoPTwoCNITraits;
#endif

private:
    // some constants from the traits for convenience
    enum {
        numEq            = BoxTraits::numEq,
        pWIdx          = TwoPTwoCNITraits::pWIdx,
        switchIdx      = TwoPTwoCNITraits::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = TwoPTwoCNITraits::temperatureIdx,
#endif

        // Phase State
        WPhaseOnly = TwoPTwoCNITraits::wPhaseOnly,
        NPhaseOnly = TwoPTwoCNITraits::nPhaseOnly,
        BothPhases = TwoPTwoCNITraits::bothPhases,

        // Grid and world dimension
        dim  = DomainTraits::dim,
        dimWorld = DomainTraits::dimWorld
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

    enum Episode {
        ImbibEpisode,  // an episode where imbibition of water takes place
        DrainEpisode,  // an episode where drainage of water takes place
        WaitEpisode    // an episode with neither drainage nor imbibition
    }; // the type of an episode of the simulation

    typedef Dune::TimeManager<Episode>                  TimeManager;
    typedef Dune::NewImplicitEulerStep<ThisType>        TimeIntegration;

    typedef typename Model::NewtonMethod                NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod>    NewtonController;

public:
    TimeDependentProblem(Grid *grid,
                         Scalar dtInitial,
                         Scalar tEnd)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          multicomp_(wPhase_, nPhase_),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("out_timeproblem")
    {
        initialTimeStepSize_ = dtInitial;
        endTime_ = tEnd;

        // specify the grid dimensions
        eps_    = 1e-8;

        depthBOR_ = 1000.0;

        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }

    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // set the episode length and initial time step size
        timeManager_.startNextEpisode(DrainEpisode, 3e2);
        timeManager_.setStepSize(initialTimeStepSize_);

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
        // execute the time integration (i.e. Runge-Kutta
        // or Euler).  TODO/FIXME: Note that the time
        // integration modifies the curTimeStepSize_ if it
        // thinks that's appropriate. (IMHO, this is an
        // incorrect abstraction, the simulation
        // controller is responsible for adapting the time
        // step size!)
        timeIntegration_.execute(*this,
                                 timeManager_.time(),
                                 stepSize,
                                 nextStepSize,
                                 1e100, // firstDt or maxDt, TODO: WTF?
                                 endTime_,
                                 1.0); // CFL factor (not relevant since we use implicit euler)
    };

    //! called by the TimeIntegration::execute function to let
    //! time pass.
    void updateModel(Scalar &dt, Scalar &nextDt)
    {
        model_.update(dt, nextDt, newtonMethod_, newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        std::cout << "Writing result file for current time step\n";

        // write the current result to disk
        writeCurrentResult_();

        // update the domain with the current solution
        //                updateDomain_();

        // change the episode of the simulation if necessary
        updateEpisode_();

        // stop the simulation if reach the end specified time
        if (timeManager_.time() >= endTime_)
            timeManager_.setFinished();
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

    void boundaryTypes(BoundaryTypeVector &values,
                       const Element &element,
                       const IntersectionIterator &isIt,
                       const GlobalPosition &globalPos,
                       const LocalPosition &localPos) const
    {
        if(globalPos[0] < eps_ || globalPos[0] > 1.0 - eps_)
            values = BoundaryConditions::neumann;
        else
            values = BoundaryConditions::dirichlet;

#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    void dirichlet(SolutionVector &values,
                   const Element &element,
                   int vertIdx,
                   int globalVertexIdx)
    {
        const LocalPosition &localPos = DomainTraits::referenceElement(element.type()).position(vertIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(vertIdx);

        // uses the same values, that are used in the initial conditions
        initial(values,
                element,
                globalPos,
                localPos);

        if(globalPos[1] > eps_) // upper boundary
        {
            if (timeManager_.episode() == DrainEpisode)
                // drain water
                values[pWIdx] = 1e5;// + (depthBOR_ - globalPos[1])*densityW*9.81;
            else if (timeManager_.episode() == ImbibEpisode) {
                // imbibition of water
                values[pWIdx] = 1.2e5;
            }
        }
    }

    /////////////////////////////
    // NEUMANN boundaries
    /////////////////////////////
    void neumann(SolutionVector &values,
                 const Element &element,
                 const IntersectionIterator &isIt,
                 const GlobalPosition &globalPos,
                 const LocalPosition &localPos) const
    {
        values = 0;
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
    void initial(SolutionVector &values,
                 const Element& element,
                 const GlobalPosition &globalPos,
                 const LocalPosition &localPos)
    {
        //                Scalar densityW = 1000.0;
        values[pWIdx] = 1e5;// + (depthBOR_ - globalPos[1])*densityW*9.81;

        if (globalPos[1] > 1.0 - eps_)
            values[switchIdx] = 0.6;
        else
            values[switchIdx] = 0.1;

#if !ISOTHERMAL
        values[temperatureIdx] = 283.15 + (depthBOR_ - globalPos[1])*0.03;
#endif
    }

    int initialPhaseState(const Vertex       &vert,
                          int              &globalIdx,
                          const GlobalPosition &globalPos) const
    {
        return BothPhases;
    }

    Scalar porosity(const Element &element, int localIdx) const
    {
        // TODO/HACK: porosity should be defined on the verts
        // as it is required on the verts!
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().porosity(globalPos, *(ParentType::elementBegin()), local);
    };

    Scalar pC(Scalar satW, int globalIdx, const GlobalPosition &globalPos)
    {
        const LocalPosition &local =
            DomainTraits::referenceElement(ParentType::elementBegin()->type()).position(0, dim);
        return materialLaw().pC(satW, globalPos, *(ParentType::elementBegin()), local);
    };


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

private:
    void updateEpisode_()
    {
        Scalar    epiLength = timeManager_.episodeLength();
        Episode   episode = timeManager_.episode();
        int       epiIdx = timeManager_.episodeIndex();

        if (timeManager_.time() >= endTime_) {
            timeManager_.setFinished();
            return;
        }

        if (!timeManager_.episodeIsOver())
            return;
        std::cerr << "episode is over!!\n";
        switch (epiIdx) {
        case 1:
            timeManager_.startNextEpisode(ImbibEpisode, 4e2);
            return;
        case 2:
            timeManager_.startNextEpisode(DrainEpisode, 3e2);
            return;
            //                    case 3:
            //                        timeManager_.startNextEpisode(WaitEpisode, 1e2);
            //                        return;
            //                    case 4:
            //                        timeManager_.startNextEpisode(DrainEpisode, 3e2);
            return;
        }

        timeManager_.startNextEpisode(episode, epiLength);
    }



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
    Scalar eps_;
    GlobalPosition  gravity_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;
    Multicomp       multicomp_;

    TimeManager     timeManager_;
    TimeIntegration timeIntegration_;
    Scalar          initialTimeStepSize_;
    Scalar          endTime_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
};
} //end namespace

#endif
