#ifndef DUNE_NIMASPROBLEM_HH
#define DUNE_NIMASPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include <dumux/material/phaseproperties/phaseproperties_waterair.hh>
#include <nimassoil.hh>
#include <dumux/material/twophaserelations.hh>

#include <dumux/material/multicomponentrelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>
#include <dune/common/collectivecommunication.hh>

#include<dumux/new_models/2p2c/2p2cboxmodel.hh>
#include<dumux/new_models/2p2cni/2p2cniboxmodel.hh>
#include<dumux/new_models/2p2cni/2p2cninewtoncontroller.hh>

#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

// with this switch, either the 2p2c or the 2p2cni model can be chosen
#define ISOTHERMAL 1

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
template<class GridT, class ScalarT>
class NimasProblem : public BasicDomain<GridT,ScalarT>
{
    typedef GridT                          			Grid;
    typedef BasicDomain<Grid, ScalarT>     			ParentType;
    typedef NimasProblem<Grid, ScalarT>    			ThisType;
    typedef TwoPTwoCPnSwTraits<ScalarT>    			Formulation;
#if !ISOTHERMAL
    typedef TwoPTwoCNIBoxModel<ThisType, Formulation>   Model;
#else
    typedef TwoPTwoCBoxModel<ThisType, Formulation>     Model;
#endif

    typedef Dune::Liq_WaterAir                     WettingPhase;
    typedef Dune::Gas_WaterAir                     NonwettingPhase;
    typedef Dune::NimasSoil<Grid, ScalarT>         Soil;
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
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

private:
    // some constants from the traits for convenience
    enum {
        numEq          = BoxTraits::numEq,
        pressureIdx    = TwoPTwoCNITraits::pressureIdx,
        switchIdx      = TwoPTwoCNITraits::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = TwoPTwoCNITraits::temperatureIdx,
#endif
        wComp = TwoPTwoCNITraits::wComp,
        nComp = TwoPTwoCNITraits::nComp,

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

    enum Episode // the type of an episode of the simulation
    {
        ConstantRateEpisode,  // an episode where imbibition of water takes place
        FallingRateEpisode,  // an episode where drainage of water takes place
        WaitEpisode    // an episode with neither drainage nor imbibition
    };

    typedef Dune::TimeManager<Episode>                  TimeManager;

    typedef typename Model::NewtonMethod                NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod>    NewtonController;

public:
    NimasProblem(Grid *grid,
                 Scalar dtInitial,
                 Scalar tEnd)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          multicomp_(wPhase_, nPhase_),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("out_nimasproblem_refined100_Xwn005")
    {
        timeManager_.setStepSize(dtInitial);
        endTime_ = tEnd;
        // specify the grid dimensions
        eps_    = 1e-8;
        depthBOR_ = 0.25;
        gravity_[0] = 0;
        //        gravity_[1] = 0; // NO gravity
        gravity_[1] = -9.81;
    }


    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

        values = BoundaryConditions::neumann;

        if (timeManager_.episode() == ConstantRateEpisode)
        {
        	if(globalPos[1] > depthBOR_ - eps_)
        	{
					values[pressureIdx] = BoundaryConditions::dirichlet;
					values[switchIdx] = BoundaryConditions::neumann;
        	}
        }
        else if (timeManager_.episode() == FallingRateEpisode)
        {
        	if(globalPos[1] > depthBOR_ - eps_)
        	{
					values[pressureIdx] = BoundaryConditions::dirichlet;
                    values[switchIdx] = BoundaryConditions::dirichlet;
        	}
        }
#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    void dirichlet(SolutionVector      &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0;
        // uses the same values, that are used in the initial conditions
        initial_(values, globalPos);

        if (timeManager_.episode() == FallingRateEpisode)
        {
        	if(globalPos[1] > depthBOR_ - eps_)
					values[switchIdx] = 0.0025;//1e-5;
//        	if(globalPos[1] > depthBOR_ - eps_)
//					values[switchIdx] = 0.032;
        }
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
            = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

        values = 0;

        if (globalPos[1] > depthBOR_ - eps_)
            {
                if (timeManager_.episode() == ConstantRateEpisode)
                    values[wComp] = 5.6e-5;
//                else if (timeManager_.episode() == FallingRateEpisode)
//                    values[wComp] = 1e-5;//5.8e-12;
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
        values = Scalar(0.0);
    }

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    void initial(SolutionVector &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

private:
	// the internal method for the initial condition
    void initial_(SolutionVector         &values,
                  const GlobalPosition   &globalPos) const
    {
        Scalar densityW = 1000.0;

        values[pressureIdx] = 1.5e5 - (depthBOR_ - globalPos[1]) * densityW * gravity_[1];
        values[switchIdx] = 1.0 - 1e-2;
        //        std::cout << globalPos[0] << ", " << globalPos[1] << ": " << values[pressureIdx] << std::endl;

#if !ISOTHERMAL
        values[temperatureIdx] = 283.15 + (depthBOR_ - globalPos[1])*0.03;
#endif
    }

public:

    int initialPhaseState(const Vertex         &vert,
                          int                  &globalIdx,
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
    };

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

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // set the initial condition
        timeManager_.startNextEpisode(ConstantRateEpisode, 822358.5);//8.05e5);//1.15e6);//1036800);
        model_.initial();

        // write the inital solution to disk
        //                writeCurrentResult_(); // TODO
    }

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
        std::cout << "episode is over!!\n";
        switch (epiIdx) {
        case 1:
            timeManager_.startNextEpisode(FallingRateEpisode, 2.5e7);//2592000);
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

        Dune::FieldVector<Scalar, 4> mass;
        model_.calculateMass(mass);
        double curTime = timeManager_.time();

        if(collectiveCom_.rank() == 0)
        {
            std::cout<< "Mass[kg] (nC, nC in nP, wC, wC in wP): "
            << curTime << ", "<< mass[0] <<", "<<mass[1]<<", "<< mass[2]<<", "<< mass[3] <<". \n";
        }

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
    Scalar          endTime_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
    CollectiveCommunication collectiveCom_;
};
} //end namespace

#endif
