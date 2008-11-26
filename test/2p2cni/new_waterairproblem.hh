#ifndef DUNE_NEW_WATERAIRPROBLEM_HH
#define DUNE_NEW_WATERAIRPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if !HAVE_UG
#error "The UG grid manager is required for this problem"
#endif


#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include <dumux/material/phaseproperties/phaseproperties_waterair.hh>
#include <dumux/material/matrixproperties.hh>
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

/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
    //! class that defines the parameters of an air waterair under a low permeable layer
    /*! Problem definition of an air injection under a low permeable layer. Air enters the domain
     * at the right boundary and migrates upwards.
     * Problem was set up using the rect2d.dgf grid.
     *
     *	Template parameters are:
     *
     *	- ScalarT  Floating point type used for scalars
     */
    template<class ScalarT>
    class NewWaterAirProblem : public BasicDomain<Dune::SGrid<2,2>, 
                                                  ScalarT>
    {
        typedef Dune::SGrid<2,2>               Grid;
        typedef BasicDomain<Grid, ScalarT>     ParentType;
        typedef NewWaterAirProblem<ScalarT>    ThisType;
        typedef TwoPTwoCNIBoxModel<ThisType>   Model;

        typedef Dune::GridPtr<Grid>                    GridPointer;

        typedef Dune::Liq_WaterAir                     WettingPhase;
        typedef Dune::Gas_WaterAir                     NonwettingPhase;
        typedef Dune::HomogeneousSoil<Grid, ScalarT>   Soil;
        typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;
        typedef Dune::CWaterAir                        Multicomp;

    public:
        // the domain traits of the domain
        typedef typename ParentType::DomainTraits   DomainTraits;
        // the traits of the BOX scheme
        typedef typename Model::BoxTraits           BoxTraits;
        // the traits of the Pw-Sn model
        typedef typename Model::TwoPTwoCNITraits    TwoPTwoCNITraits;

    private:
        // some constants from the traits for convenience
        enum {
            PrimaryVariables = BoxTraits::PrimaryVariables,
            PwIndex          = TwoPTwoCNITraits::PwIndex,
            SwitchIndex      = TwoPTwoCNITraits::SwitchIndex,
            TemperatureIndex = TwoPTwoCNITraits::TemperatureIndex,
            
            // Phase State
            WPhaseOnly = TwoPTwoCNITraits::WPhaseOnly,
            NPhaseOnly = TwoPTwoCNITraits::NPhaseOnly,
            BothPhases = TwoPTwoCNITraits::BothPhases,

            // Grid and world dimension
            GridDim  = DomainTraits::GridDim,
            WorldDim = DomainTraits::WorldDim
        };
      
        // copy some types from the traits for convenience
        typedef typename DomainTraits::Scalar                     Scalar;
        typedef typename DomainTraits::Cell                       Cell;
        typedef typename DomainTraits::CellIterator               CellIterator;
        typedef typename DomainTraits::ReferenceElement           ReferenceElement;
        typedef typename DomainTraits::Node                       Node;
        typedef typename DomainTraits::NodeIterator               NodeIterator;
        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter IntersectionIteratorGetter;
        typedef typename DomainTraits::LocalCoord                 LocalCoord;
        typedef typename DomainTraits::WorldCoord                 WorldCoord;

        typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction               SpatialFunction;
        typedef typename BoxTraits::UnknownsVector                UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

        typedef VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

        enum Episode {}; // the type of an episode of the simulation
        typedef TimeManager<Episode>                        TimeManager;
        typedef Dune::NewImplicitEulerStep<ThisType>        TimeIntegration;

        typedef typename Model::NewtonMethod                NewtonMethod;
        typedef TwoPTwoCNINewtonController<NewtonMethod>    NewtonController;

    public:
        NewWaterAirProblem(Grid *grid,
                           Scalar dtInitial,
                           Scalar tEnd) 
            : ParentType(grid),
              materialLaw_(soil_, wPhase_, nPhase_),
              multicomp_(wPhase_, nPhase_),
              model_(*this),
              newtonMethod_(model_),
              resultWriter_("new_waterair")
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
                timeManager_.startNextEpisode(1e100);
                timeManager_.setStepSize(initialTimeStepSize_);

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
                writeCurrentResult_(); // TODO

                // update the domain with the current solution
//                updateDomain_();

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
          \return	wetting phase
        */
        const WettingPhase &wettingPhase() const
            { return wPhase_; }

        //! properties of the nonwetting (liquid) phase
        /*! properties of the nonwetting (liquid) phase
          \return	nonwetting phase
        */
        const NonwettingPhase &nonwettingPhase() const
            { return nPhase_; }
 
          
        //! properties of the soil
        /*! properties of the soil
          \return	soil
        */
        const Soil &soil() const
            {  return soil_; }

        //! properties of the soil
        /*! properties of the soil
          \return	soil
        */
        Soil &soil()
            {  return soil_; }

        //! object for multicomponent calculations
        /*! object for multicomponent calculations including mass fractions,
         * mole fractions and some basic laws
         \return	multicomponent object
        */
        MultiComp &multicomp ()
//        const MultiComp &multicomp () const
            {
                return multicomp_;
            }

        //! object for definition of material law
        /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
          \return	material law
        */
        MaterialLaw &materialLaw ()
//        const MaterialLaw &materialLaw () const
            {
                return materialLaw_;
            }
        
        void boundaryTypes(BoundaryTypeVector &values,
                           const Cell &cell,
                           const IntersectionIterator &faceIt,
                           const WorldCoord &globalPos,
                           const LocalCoord &localPos) const
            {
                if(globalPos[0] < eps_)
                    values = BoundaryConditions::dirichlet;
                else
                    values = BoundaryConditions::neumann;
            }

        /////////////////////////////
        // DIRICHLET boundaries
        /////////////////////////////
        void dirichlet(UnknownsVector &values,
                       const Cell &cell,
                       int nodeIdx,
                       int globalNodeIdx)
            {
                const LocalCoord &localPos = DomainTraits::referenceElement(cell.type()).position(nodeIdx, GridDim);
                const WorldCoord &globalPos = cell.geometry()[nodeIdx];
                
                initial(values,
                        cell,
                        globalPos, 
                        localPos);
            }

        /////////////////////////////
        // NEUMANN boundaries
        /////////////////////////////
        void neumann(UnknownsVector &values,
                     const Cell &cell,
                     const IntersectionIterator &faceIt,
                     const WorldCoord &globalPos,
                     const LocalCoord &localPos) const
            {
                values = 0;

                // negative values for injection
                if (globalPos[1] > 1.0 && globalPos[1] < 5.0)
                {
                    values[SwitchIndex] = -1e-5;
                }
            }

        /////////////////////////////
        // sources and sinks
        /////////////////////////////
        void sourceTerm(UnknownsVector &values,
                        const Cell &cell,
                        const FVElementGeometry &, 
                        int subControlVolumeIdx) const
            {
                values = Scalar(0.0);
            }

        //////////////////////////////
      
        /////////////////////////////
        // INITIAL values
        /////////////////////////////
        void initial(UnknownsVector &values,
                     const Cell& cell,
                     const WorldCoord &globalPos, 
                     const LocalCoord &localPos)
            {
                Scalar densityW = 1000.0;
                values[PwIndex] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
                values[SwitchIndex] = 0.0;
                values[TemperatureIndex] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
            }

        Scalar porosity(const Cell &cell, int localIdx) const
            {
                // TODO/HACK: porosity should be defined on the nodes
                // as it is required on the nodes!
                const LocalCoord &local =
                    DomainTraits::referenceElement(cell.type()).position(localIdx, GridDim);
                const WorldCoord &globalPos = cell.geometry()[localIdx];
                return soil().porosity(globalPos, *(ParentType::cellBegin()), local);
            };

        Scalar pC(Scalar satW, int globalIdx, const WorldCoord &globalPos)
            {
                // TODO/HACK: porosity should be defined on the nodes
                // as it is required on the nodes!
                const LocalCoord &local =
                    DomainTraits::referenceElement(ParentType::cellBegin()->type()).position(0, GridDim);
                return materialLaw().pC(satW, globalPos, *(ParentType::cellBegin()), local);
            };


        int initialPhaseState(const Node       &node,
                              int              &globalIdx,
                              const WorldCoord &globalPos) const
            {
                return WPhaseOnly;
            }

      
        const WorldCoord &gravity () const
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
        WorldCoord  gravity_;

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
