#ifndef DUNE_NEW_INJECTIONPROBLEM_HH
#define DUNE_NEW_INJECTIONPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
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

#include<dumux/new_models/2p2c/2p2cboxmodel.hh>
#include<dumux/new_models/2p2c/2p2cnewtoncontroller.hh>

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
//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
    template<class Grid, class ScalarT>
    class InjectionSoil: public Matrix2p<Grid,ScalarT>
    {
    public:
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef ScalarT Scalar;
        typedef typename Grid::ctype CoordScalar;
        enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

        typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
        typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

        InjectionSoil():Matrix2p<Grid,Scalar>()
            {
                lowK_ = highK_ = 0.;
                for(int i = 0; i < dim; i++){
                    lowK_[i][i] = 1e-13;
                    highK_[i][i] = 1e-12;
                }
                layerBottom_ = 22.0;
            }

        ~InjectionSoil()
            {}

        const FieldMatrix<CoordScalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi)
            {
                if (x[1] < layerBottom_)
                    return highK_;
                else
                    return lowK_;
            }

        double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
            {
                return 0.3;
            }

        double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
            {
                return 0.2;
            }

        double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
            {
                return 0.0;
            }

        /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
         * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
        double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
            {
                return     790 /* spec. heat cap. of granite */
                    * 2700 /* density of granite */
                    * porosity(x, e, xi);
            }

        double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
            {
                static const double lWater = 0.6;
                static const double lGranite = 2.8;
                double poro = porosity(x, e, xi);
                double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
                double ldry = pow(lGranite, (1-poro));
                return ldry + sqrt(sat) * (ldry - lsat);
            }

        std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
            {
                // example for Brooks-Corey parameters
                std::vector<double> param(2);
                param[0] = 2.; // lambda
                param[1] = 1e4; // entry-pressures

                return param;
            }

        typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
            {
                return Matrix2p<Grid,Scalar>::brooks_corey;
            }

    private:
        FieldMatrix<Scalar,dim,dim> lowK_, highK_;
        double layerBottom_;
    };

    //! class that defines the parameters of an air injection under a low permeable layer
    /*! Problem definition of an air injection under a low permeable layer. Air enters the domain
     * at the right boundary and migrates upwards.
     * Problem was set up using the rect2d.dgf grid.
     *
     *    Template parameters are:
     *
     *    - ScalarT  Floating point type used for scalars
     */
    template<class ScalarT>
    class NewInjectionProblem : public BasicDomain<Dune::YaspGrid<2>,
                                                   ScalarT>
    {
        typedef Dune::YaspGrid<2>              Grid;
        typedef BasicDomain<Grid, ScalarT>     ParentType;
        typedef NewInjectionProblem<ScalarT>   ThisType;
        typedef TwoPTwoCBoxModel<ThisType>     Model;

        typedef Dune::GridPtr<Grid>                    GridPointer;

        typedef Dune::Liq_WaterAir                     WettingPhase;
        typedef Dune::Gas_WaterAir                     NonwettingPhase;
        typedef Dune::InjectionSoil<Grid, ScalarT>     Soil;
        typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;
        typedef Dune::CWaterAir                        Multicomp;

    public:
        // the domain traits of the domain
        typedef typename ParentType::DomainTraits   DomainTraits;
        // the traits of the BOX scheme
        typedef typename Model::BoxTraits           BoxTraits;
        // the traits of the Pw-Sn model
        typedef typename Model::TwoPTwoCTraits      TwoPTwoCTraits;

    private:
        // some constants from the traits for convenience
        enum {
            numEq       = BoxTraits::numEq,
            pWIdx     = TwoPTwoCTraits::pWIdx,
            switchIdx = TwoPTwoCTraits::switchIdx,

            // Phase State
            WPhaseOnly = TwoPTwoCTraits::wPhaseOnly,
            NPhaseOnly = TwoPTwoCTraits::nPhaseOnly,
            BothPhases = TwoPTwoCTraits::bothPhases,

            // Grid and world dimension

            dim      = DomainTraits::dim,
            dimWorld = DomainTraits::dimWorld
        };

        // copy some types from the traits for convenience
        typedef typename DomainTraits::Scalar                     Scalar;
        typedef typename DomainTraits::Element                       Element;
        typedef typename DomainTraits::ElementIterator               ElementIterator;
        typedef typename DomainTraits::ReferenceElement           ReferenceElement;
        typedef typename DomainTraits::Vertex                       Vertex;
        typedef typename DomainTraits::VertexIterator               VertexIterator;
        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::LocalPosition                 LocalPosition;
        typedef typename DomainTraits::GlobalPosition                 GlobalPosition;

        typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction               SpatialFunction;
        typedef typename BoxTraits::SolutionVector                SolutionVector;
        typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

        typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

        enum Episode {}; // the type of an episode of the simulation
        typedef Dune::TimeManager<Episode>                  TimeManager;
        typedef Dune::NewImplicitEulerStep<ThisType>        TimeIntegration;

        typedef typename Model::NewtonMethod                NewtonMethod;
        typedef TwoPTwoCNewtonController<NewtonMethod>      NewtonController;

    public:
        NewInjectionProblem(Grid *grid,
                            Scalar dtInitial,
                            Scalar tEnd)
            : ParentType(grid),
              materialLaw_(soil_, wPhase_, nPhase_),
              multicomp_(wPhase_, nPhase_),
              timeManager_(this->grid().comm().rank() == 0),
              model_(*this),
              newtonMethod_(model_),
              resultWriter_("new2p2c")
            {
                initialTimeStepSize_ = dtInitial;
                endTime_ = tEnd;

                // specify the grid dimensions
                outerLowerLeft_[0] = 0;
                outerLowerLeft_[1] = 0;
                outerUpperRight_[0] = 50.0;
                outerUpperRight_[1] = 40.0;

                height_ = outerUpperRight_[1] - outerLowerLeft_[1];
                width_  = outerUpperRight_[0] - outerLowerLeft_[0];
                eps_    = 1e-8*width_;

                // for defining e.g. a lense
                innerLowerLeft_[0] = 0.0;
                innerLowerLeft_[1] = 0.0;

                innerUpperRight_[0] = 0.0;
                innerUpperRight_[1] = 0.5;

                depthBOR_ = 800.0;

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
                if (this->grid().comm().rank() == 0)
                    std::cout << "Writing result file for current time step\n";

                // write the current result to disk
                writeCurrentResult_();

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
                if (globalPos[0] < eps_)
                    values = BoundaryConditions::dirichlet;
                else
                    values = BoundaryConditions::neumann;

//        if (globalPos[1] < eps_)
//            values = BoundaryConditions::dirichlet;
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

                initial(values,
                        element,
                        globalPos,
                        localPos);
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

                //Scalar lambda = (globalPos[1])/height_;
                if (globalPos[1] < 15 && globalPos[1] > 5) {
                    values[switchIdx] = -1e-3;
                }
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
                Scalar densityW_ = 1000.0;

                values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - globalPos[1]);
                values[switchIdx] = 1e-8;

//                std::cout << "element " << ParentType::elementIdx(element) << " position of vert " << globalVertexIdx << ": " << globalPos << " -> " << values << "\n";

//        if ((globalPos[0] > 60.0 - eps_) && (globalPos[1] < 10 && globalPos[1] > 5))
//            values[switchIdx] = 0.05;

//            if (globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]
//             && globalPos[0] >= innerLowerLeft_[0])
//                values[switchIdx] = 0.2;
//            else
//                values[switchIdx] = 1e-6;
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
                // TODO/HACK: porosity should be defined on the verts
                // as it is required on the verts!
                const LocalPosition &local =
                    DomainTraits::referenceElement(ParentType::elementBegin()->type()).position(0, dim);
                return materialLaw().pC(satW, globalPos, *(ParentType::elementBegin()), local);
            };


        int initialPhaseState(const Vertex       &vert,
                              int              &globalIdx,
                              const GlobalPosition &globalPos) const
            {
                int state;

//        state = BothPhases;
                state = WPhaseOnly;

//        if ((globalPos[0] > 60.0 - eps_) && (globalPos[1] < 10 && globalPos[1] > 5))
//            state = bothPhases;
//            if (globalPos[1] >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]
//                  && globalPos[0] >= innerLowerLeft_[0])
//                state = 2;
//            else

                return state;
            }


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
        Scalar depthBOR_;
        Scalar eps_;
//      Scalar densityW_, densityN_;
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
