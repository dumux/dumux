#ifndef DUNE_NEW_BLOBPROBLEM_HH
#define DUNE_NEW_BLOBPROBLEM_HH

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

#include<dumux/new_models/2p2c/2p2cboxmodel.hh>
#include<dumux/new_models/2p2c/2p2cnewtoncontroller.hh>

#include<dumux/timedisc/new_impliciteulerstep.hh>
#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

/**
 * @file
 *
 * @brief Definition of a problem, where a blob of gas is enclosed by
 * a zone completely saturated with water. The gas saturation within
 * the blob is below the residual saturation, bit gas gets transported
 * away anyway because it is partially miscible with water.
 *
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
    template<class Grid, class ScalarT>
    class BlobSoil: public Matrix2p<Grid,ScalarT>
    {
    public:
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef ScalarT Scalar;

        enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

        typedef Dune::FieldVector<Scalar,dim>      LocalPosition;
        typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

        BlobSoil()
            : Matrix2p<Grid,Scalar>(),
              K_(0.0)
            {
                for(int i = 0; i < dim; i++)
                    K_[i][i] = 1e-12;
            }

        ~BlobSoil()
            {}

      virtual const FieldMatrix<Scalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi)
      {
          return K_;
      }
      virtual double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
      {
          return 0.3;
      }

      virtual double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
      {
          return 0;
      }

      virtual double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
      {
          return 0.1;
      }

      /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
               * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
      virtual double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
      {
          return     790 /* spec. heat cap. of granite */
                          * 2700 /* density of granite */
                          * (1 - porosity(x, e, xi));
      }

      virtual double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
      {
          static const double lWater = 0.6;
          static const double lGranite = 2.8;
          double poro = porosity(x, e, xi);
          double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
          double ldry = pow(lGranite, (1-poro));
          return ldry + sqrt(sat) * (ldry - lsat);
      }

      virtual std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
      {
          // example for Brooks-Corey parameters
          std::vector<double> param(2);
          param[0] = 2.; // lambda
          param[1] = 0.; // entry-pressures

  //        if (x[0] > 150)
  //            param[0] = 0.5;

          return param;
      }

      virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
      {
          return Matrix2p<Grid,Scalar>::brooks_corey;
      }

  private:
      FieldMatrix<Scalar,dim,dim> K_;
    };

    /*!
     * @brief Definition of a problem, where a blob of gas is enclosed by
     * a zone completely saturated with water. The gas saturation within
     * the blob is below the residual saturation, bit gas gets transported
     * away anyway because it is partially miscible with water.
     *
     *    Template parameters are:
     *
     *    - ScalarT  Floating point type used for scalars
     */
    template<class ScalarT>
    class NewBlobProblem : public BasicDomain<Dune::SGrid<2, 2>,
                                               ScalarT>
    {
        typedef Dune::SGrid<2,2>               Grid;
        typedef BasicDomain<Grid, ScalarT>     ParentType;
        typedef NewBlobProblem<ScalarT>        ThisType;
        typedef TwoPTwoCBoxModel<ThisType>     Model;

        typedef Dune::GridPtr<Grid>                    GridPointer;

        typedef Dune::Liq_WaterAir                     WettingPhase;
        typedef Dune::Gas_WaterAir                     NonwettingPhase;
        typedef Dune::BlobSoil<Grid, ScalarT>          Soil;
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
            wPhaseOnly = TwoPTwoCTraits::wPhaseOnly,
            nPhaseOnly = TwoPTwoCTraits::nPhaseOnly,
            bothPhases = TwoPTwoCTraits::bothPhases,

            // Grid and world dimension
            dim  = DomainTraits::dim,
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
        NewBlobProblem(Scalar dtInitial,
                       Scalar tEnd)
            : ParentType(new Grid(Dune::FieldVector<int,dim>(40), // number of verts
                                  GlobalPosition(0.0),  // lower left
                                  GlobalPosition(300.0) // upper right
                             )),
              materialLaw_(soil_, wPhase_, nPhase_),
              multicomp_(wPhase_, nPhase_),
              model_(*this),
              newtonMethod_(model_),
              resultWriter_("newblob")
            {
                initialTimeStepSize_ = dtInitial;
                endTime_ = tEnd;

                eps_    = 1e-8 * 300.0;

                depthBOR_ = 800.0;

                gravity_[0] = 0;
                gravity_[1] = 0; // -9.81;
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
                values = BoundaryConditions::neumann;

                if ((globalPos[0] < eps_) || (globalPos[0] > (300 - eps_)))
                    values = BoundaryConditions::dirichlet;
            }

        /////////////////////////////
        // DIRICHLET boundaries
        /////////////////////////////
        void dirichlet(SolutionVector &values,
                       const Element &element,
                       int vertIdx,
                       int globalVertexIdx)
            {
//                const LocalPosition &localPos = DomainTraits::referenceElement(element.type()).position(vertIdx, dim);
                const GlobalPosition &globalPos = element.geometry().corner(vertIdx);

                values[pWIdx] = 1e5;
                values[switchIdx] = 0.0;

                if (globalPos[0] < eps_)
                {
                    values[pWIdx] = 1e5 + 50; // used to be 2e5, but then diffusion is negligible
                    values[switchIdx] = 0;  // may be Sn, Xaw or Xwn!!
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
                values = 0.0;
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
                values[pWIdx] = 1e5;//(600-globalPos[0])/300 * 1e5;
                values[switchIdx] = 0;

                if (isInsideBlob_(globalPos))
                {
                    values[switchIdx] = 0.1;
                }
            }

        int initialPhaseState(const Vertex       &vert,
                              int              &globalIdx,
                              const GlobalPosition &globalPos) const
            {
                if (isInsideBlob_(globalPos))
                    return bothPhases;

                return wPhaseOnly;
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
        bool isInsideBlob_(const GlobalPosition &globalPos) const
            {
                return globalPos[0] >= 59.0 &&
                    globalPos[0] <= 121 &&
                    globalPos[1] >= 119 &&
                    globalPos[1] <= 181;
            }


        // write the fields current solution into an VTK output file.
        void writeCurrentResult_()
            {
                resultWriter_.beginTimestep(timeManager_.time(),
                                            ParentType::grid().leafView());

                model_.addVtkFields(resultWriter_);

                resultWriter_.endTimestep();
            }


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
