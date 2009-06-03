#ifndef DUNE_SEQCOUP2P2CNI_HH
#define DUNE_SEQCOUP2P2CNI_HH

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
#include <dumux/material/phaseproperties/phaseproperties_brineco2.hh>

#include <dumux/material/multicomponentrelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/io.hh>

#include<dumux/new_models/2p2cni/2p2cniboxmodel.hh>
#include<dumux/new_models/2p2c/2p2cboxjacobianbase.hh>
#include <dumux/new_models/2p2c/2p2cvertexdata.hh>
#include<dumux/new_models/2p2cni/2p2cninewtoncontroller.hh>


#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>
#include <dumux/io/restart.hh>

#define ISOTHERMAL 0

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
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template<class ScalarT>
class SeqCoup2P2CNIProblem : public BasicDomain<Dune::YaspGrid<2>,
                                              ScalarT>
{
    typedef Dune::YaspGrid<2>              Grid;
    typedef BasicDomain<Grid, ScalarT>     ParentType;
    typedef SeqCoup2P2CNIProblem<ScalarT>    ThisType;
#if !ISOTHERMAL
    typedef TwoPTwoCNIBoxModel<ThisType>   Model;
#else
    typedef TwoPTwoCBoxModel<ThisType>   Model;
#endif

    typedef Dune::GridPtr<Grid>                    GridPointer;

    typedef Dune::Liq_BrineCO2                     WettingPhase;
    typedef Dune::Gas_BrineCO2                     NonwettingPhase;
    typedef Dune::HomogeneousSoil<Grid, ScalarT>   Soil;
    typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;
    typedef Dune::CBrineCO2                        Multicomp;

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
    typedef typename Dune::BlockVector<Dune::FieldVector<double,3> > SolVector;
    typedef typename Dune::BlockVector<Dune::FieldVector<double,1> > PhaseState;
    typedef typename Dune::BlockVector<Dune::FieldVector<bool,1> > BoolVector;
private:
    // some constants from the traits for convenience
    enum {
        numEq          = BoxTraits::numEq,
        pressureIdx    = TwoPTwoCNITraits::pressureIdx,
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
        dimWorld = DomainTraits::dimWorld,
        wPhase = 0,
        nPhase = 1,
        wComp = 0,
        nComp = 1,
        pWsN        = 0,
        pNsW        = 1
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
    typedef typename BoxTraits::LocalFunction                 LocalFunction;

    typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>                  TimeManager;

    typedef typename Model::NewtonMethod                NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod>    NewtonController;

public:
    SeqCoup2P2CNIProblem(Grid *grid,
                       Scalar dtInitial,
                       Scalar tEnd,
                       bool sequentialCoupling = true)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          multicomp_(wPhase_, nPhase_),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("seqcoup_2p2cni"),
          phaseState_(this->numVertices()),
          corrected_(this->numVertices())

    {
        timeManager_.setStepSize(dtInitial);

        // specify the grid dimensions
        eps_    = 1e-8;

        depthBOR_ = 3238.2;

        gravity_[0] = 0;
        gravity_[dim - 1] = -9.81;

        // choose primary variables
        formulation_ = pWsN;

        wasRestarted_ = false;
        sequentialCoupling_ = sequentialCoupling;
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

        // write restart file after every five steps
        static int dummy = 0;
        ++dummy;
        if (dummy % 1 == 0)
            serialize();
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
        values = BoundaryConditions::neumann;

        if (globalPos[0] > 99.)
            values = BoundaryConditions::dirichlet;

#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
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

        values = 0;
        // negative values for injection
        if (globalPos[1] < 10.0 && globalPos[0] < 1.0)
        	values[switchIdx] = -0.005;
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
        int globalIdx = this->vertexIdx(element, scvIdx);
        values = initialSolution_[globalIdx];
     }

    void sequentialCoupling(SolVector& twoPNISolution)
    {
        setPhaseState_(twoPNISolution);
        correctSaturation_(twoPNISolution);
        initialSolution_ = twoPNISolution;
    }
private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(SolutionVector       &values,
                  const GlobalPosition &globalPos) const
    {
           Scalar densityB = 1037.24;
           Scalar pRef = 101300.;
           Scalar TRef = 283.15;

            values[pressureIdx] = pRef - (depthBOR_ - globalPos[dim - 1]) * densityB*gravity_[1];
            values[switchIdx] = 0.0;
            values[temperatureIdx] = TRef + (depthBOR_ - globalPos[dim - 1]) * 0.03;
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
                           "seqcoup_2p2cni",
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

        res.serializeEnd();
    }
public:
    Scalar porosity(const Element &element, int localIdx) const
    {
        // TODO/HACK: porosity should be defined on the verts
        // as it is required on the verts!
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().porosity(globalPos, element, local);
    };

    FieldMatrix<Scalar,dim,dim> &K(const Element &element, int localIdx)
    {
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().K(globalPos, element, local);
    }

    void deserialize(double t)
    {
        typedef Dune::Restart<Grid> Restarter;
        Restarter res;
        res.deserializeBegin(this->grid(), "seqcoup_2pni", t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };


    int initialPhaseState(const Vertex       &vert,
                          int              &globalIdx,
                          const GlobalPosition &globalPos) const
    {
       return phaseState_[globalIdx];
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

    int formulation () const
    {
        return formulation_;
    }
private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        Dune::FieldVector<Scalar, 4> mass;
        std::cout<<"\nMINIMUM AND MAXIMUM VALUES:"<<std::endl;
        model_.calculateMass(mass);
        std::cout<<"\nMASS BALANCE AND TIMESTEP INFORMATIONS:"<<std::endl;
        std::cout<< "Mass CO2: "<< mass[0]<<"Mass CO2 in Phase: "<<mass[1] <<"Mass Brine: "<< mass[2]<<"Mass Brine in Phase: "
        << mass[3] <<".       ";
        resultWriter_.endTimestep();
    }

    void setPhaseState_(SolVector& twoPNISolution)
        {
            int numVertices = this->numVertices();
            for (int globalIdx = 0; globalIdx < numVertices; ++globalIdx)
            {
                if(twoPNISolution[globalIdx][switchIdx] <= 1.e-10)
                {
                    // initialize variable phaseState
                    phaseState_[globalIdx] = WPhaseOnly;
                }

                else if(twoPNISolution[globalIdx][switchIdx] >= 1-1.e-10)
                {
                    // initialize variable phaseState
                    phaseState_[globalIdx] = NPhaseOnly;
                }

                else if(twoPNISolution[globalIdx][switchIdx] < 1-1.e-10
                        || twoPNISolution[globalIdx][switchIdx]> 1.e-10)
                {
                    // initialize variable phaseState
                    phaseState_[globalIdx] = BothPhases;
                }
            }
        };

    void correctSaturation_(SolVector& twoPNISolution)
        {
            ElementIterator elementIt = this->elementBegin();
            ElementIterator endit = this->elementEnd();
            Scalar result, corrTolerance;
            Scalar xCorrNew, xCorrOld;
            Scalar satCorrNew, satCorrOld;
            Scalar twoPNIDens, twoPNISatN;
            corrTolerance = 1.e-9;
            Scalar satW, satN, pC, pN, pW, temperature;
            FieldVector<Scalar,2> density;
            FieldMatrix<Scalar,2,2> massfrac;
            corrected_ = false;
            // Loop over elements
            for (; elementIt != endit; ++elementIt)
            {
                int numLocalVerts = elementIt->template count<dim>();
                const Element& element = *elementIt;
                unsigned numVertices = this->numVertices();
                // Loop over element vertices
                for (int idx = 0; idx < numLocalVerts; ++idx)
                {
                    const LocalPosition &local = DomainTraits::referenceElement(element.type()).position(idx, dim);
                    const GlobalPosition &globalPos = element.geometry().corner(idx);
                    int globalIdx = this->vertexIdx(element, idx);

                    if(phaseState_[globalIdx] == BothPhases && corrected_[globalIdx] == false)
                    {
                        satN = twoPNISolution[globalIdx][switchIdx];
                        satW = 1 - satN;
                        temperature = twoPNISolution[globalIdx][temperatureIdx];
                        pW = twoPNISolution[globalIdx][pressureIdx];
                        pC = materialLaw().pC(satW, globalPos, element, local);
                        pN = pW + pC;
                        density[wPhase] = wPhase_.density(temperature, pW);
                        density[nPhase] = nPhase_.density(temperature, pN);
                        massfrac[nComp][wPhase] = multicomp().xAW(pN, temperature);
                        massfrac[wComp][nPhase] = multicomp().xWN(pN, temperature);
                        massfrac[wComp][wPhase] = 1.0 - massfrac[nComp][wPhase];
                        massfrac[nComp][nPhase] = 1.0 - massfrac[wComp][nPhase];
                        twoPNIDens = density[nPhase];
                        twoPNISatN = satN;
                        satCorrNew = 1.;
                        satCorrOld = 0.;

                        while(std::fabs(satCorrNew - satCorrOld)> corrTolerance && satCorrNew> 0.)// iteration loop for correct saturation value to take
                        // capillary pressure effects into account

                        {
                            satCorrOld = satCorrNew;
                            satCorrNew = (twoPNIDens * twoPNISatN
                                    - satW * density[wPhase] * massfrac[nComp][wPhase])/
                                   (density[nPhase] * massfrac[nComp][nPhase]);

                            // correction of densities and massfractions due to pressure changes
                            satW = 1.0 - satCorrNew;
                            pC = materialLaw().pC(satW, globalPos, element, local);
                            pN = pW + pC;
                            massfrac[nComp][wPhase] = multicomp().xAW(pN, temperature);
                            massfrac[wComp][nPhase] = multicomp().xWN(pN, temperature);
                            massfrac[wComp][wPhase] = 1.0 - massfrac[nComp][wPhase];
                            massfrac[nComp][nPhase] = 1.0 - massfrac[wComp][nPhase];
                            density[wPhase] = wPhase_.density(temperature, pW, massfrac[nComp][wPhase]);
                            density[nPhase] = nPhase_.density(temperature, pN);
                        }
                        result = satCorrNew;

                        // if the resulting saturation becomes < 0 all co2 of the nPhase can
                        // be dissolved in the brine. The new mass fraction is calculated here
                        if(satCorrNew <= 0.0)
                        {
                            phaseState_[globalIdx] = WPhaseOnly;

                            xCorrNew = 1.;
                            xCorrOld = 0.;
                            while(std::fabs(xCorrNew - xCorrOld)> corrTolerance)
                            {
                                xCorrOld = xCorrNew;
                                xCorrNew = (twoPNIDens * twoPNISatN)/density[wPhase];
                                density[wPhase] = wPhase_.density(temperature, pW, xCorrNew);
                            }
                            result = xCorrNew;
                        }
                        twoPNISolution[globalIdx][switchIdx] = result;
//                        std::cout<<"XCORR   "<<"XCO2n = "<<result <<"   "<<(result*density[wPhase])<<"   "<<(twoPNIDens*twoPNISatN)<<std::endl;
//                        std::cout<<"SatCORR "<<"Sn = "<<satCorrNew<<"    "<<(result*density[nPhase]*massfrac[nComp][nPhase]+(1-result)*massfrac[nComp][wPhase]*density[wPhase])<<" "<<(twoPNIDens*twoPNISatN)<<std::endl;
//                        std::cout<<"newSwitchIdx: "<<result<<", "<<twoPNISolution[globalIdx]<<std::endl;
                        corrected_[globalIdx] = true;
                    }
                }
            }
        };

     Scalar depthBOR_;
    Scalar eps_;
    int formulation_;
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
    SolVector       initialSolution_;
    PhaseState      phaseState_;
    BoolVector      corrected_;

    bool wasRestarted_;
    bool sequentialCoupling_;
};
} //end namespace

#endif
